#!/usr/bin/env Rscript
# =============================================================================
# gt_score_rediportal_v9.R   GT_REDI_PART=universe|score
#
# REDIportal v2 as truth under the v9 rules. As gt_score_rediportal_v8.R except:
#   1. universe from both read orientations: merged NTC counts at every 12-type candidate position
#      (gt_eval_v6/spectrum12: candidates per (site, strand, alt) at cov >= 5, alt >= 1, >= 1%; dbSNP common removed;
#      BQ/MQ >= 30). At a genomic A (T) the forward orientation is read A>G on + (-) reads, the reverse orientation
#      read T>C on - (+) reads. Per rule, an orientation qualifies if it passes the rule (cov >= 10, alt, edit);
#      where both qualify the one with more reads is kept (ties: forward). Multiallelic check on the kept orientation.
#      dbSNP common and GT_REDI_EXCLUDE_BED (PGP1 PASS + VQLOW) positions removed (as v8).
#   2. germline negatives: PGP1 PASS SNV negatives of the v9 KO truth (neg_source pgp1_ag / pgp1_hom,
#      gt_eval_v9/truth/gt_v9_truth_sites.tsv.gz) are negatives for every rule, in addition to U_v \ REDIportal.
#   3. per-unit recall denominator: forward-orientation positives use the REDIportal unit pileups (as v8),
#      reverse-orientation positives the reverse-strand unit pileup (GT_REDI_REVCOV/ko_truth/<bam>.tsv.gz,
#      gt_v8_unit_cov.sh LIB=reverse at rev_positive_sites).
#   call set  = AG rows (RNA frame) of the call files; a call at a reverse-orientation positive is AG only if the
#               method placed it on the A-to-I strand
#   precision = calls at positives_v / calls at (U_v + germline negatives)
#   recall, F1, tier summaries as v8
#   universe  -> GT_REDI_UDIR/{universe.rds, positives.tsv, orientation.tsv, rev_positive_sites.tsv.gz,
#                positive_sites_<v>.tsv.gz, removed_sites.tsv}
#   score     -> GT_REDI_OUT/{unit_scores.tsv, gt_tier_scoring_<v>.tsv} (+ links to the universe tables)
# Env as v8 (GT_RES_DIRS, GT_CALL_RE, GT_UNIT_MIN_COV, GT_UNIT_KEEP, GT_REDI_EXCLUDE_BED) + GT_REDI_UDIR, GT_REDI_REVCOV
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("SCRNA_DIR", ".")
PM <- file.path(B, "human/pm")
V6 <- file.path(PM, "gt_eval_v6"); V9 <- file.path(PM, "gt_eval_v9")
RD <- file.path(V6, "rediportal"); SP <- file.path(V6, "spectrum12")
PART <- Sys.getenv("GT_REDI_PART", "score")
UDIR <- Sys.getenv("GT_REDI_UDIR", file.path(V9, "rediportal_pwcov", "universe"))
REVCOV <- Sys.getenv("GT_REDI_REVCOV", file.path(V9, "rediportal_pwcov", "unit_cov_rev"))
CHRS <- c(paste0("chr", 1:22), "chrX", "chrY")
RULES <- list(cov10_alt1 = c(alt = 1, edit = 0.01), alt3 = c(alt = 3, edit = 0.01), alt3_e05 = c(alt = 3, edit = 0.05))
key <- function(chr, pos) match(chr, CHRS) * 1e9 + as.numeric(pos)
site_key <- function(s) {
  if (!length(s)) return(numeric())
  p <- tstrsplit(s, ":", fixed = TRUE); key(p[[1]], p[[2]])
}
key_site <- function(k) sprintf("%s:%d", CHRS[k %/% 1e9], as.integer(k %% 1e9))
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

if (PART == "universe") {
  dir.create(UDIR, recursive = TRUE, showWarnings = FALSE)
  redi <- rbindlist(lapply(CHRS, function(c)
    fread(file.path(RD, "bed", paste0(c, ".bed")), header = FALSE, select = c(1, 3), col.names = c("chr", "pos"))))
  redi_key <- unique(key(redi$chr, redi$pos)); rm(redi)
  snp_key <- unique(c(site_key(fread(file.path(RD, "rediportal_dbsnp_common_sites.txt"), header = FALSE)$V1),
                      site_key(fread(file.path(V6, "candidates", "candidates_annot.tsv"))[in_dbsnp_common == 1, site])))

  fs <- list.files(file.path(SP, "pileups"), "^merged_NTC_ctrl__.*\\.tsv\\.gz$", full.names = TRUE)
  stopifnot(length(fs) == 90L)
  p <- rbindlist(lapply(fs, function(f) fread(f, showProgress = FALSE)))
  p[, ref := toupper(ref)]
  p <- p[ref %in% c("A", "T") & chr %in% CHRS]
  p[, k := key(chr, pos)]
  p <- unique(p, by = "k")
  isA <- p$ref == "A"
  # forward = A-to-I strand (A: + reads, T: - reads); reverse = the opposite strand
  p[, `:=`(alt_f = fifelse(isA, G, c), cov_f = fifelse(isA, ref_plus + G, ref_minus + c),
           t1_f = fifelse(isA, C, a), t2_f = fifelse(isA, T, g),
           alt_r = fifelse(isA, g, C), cov_r = fifelse(isA, ref_minus + g, ref_plus + C),
           t1_r = fifelse(isA, c, A), t2_r = fifelse(isA, t, G))]
  p <- p[alt_f > 0 | alt_r > 0, .(k, alt_f, cov_f, t1_f, t2_f, alt_r, cov_r, t1_r, t2_r)]
  rm(isA)
  msg("merged NTC A/T candidate positions: ", nrow(p))

  drop_log <- list(dbsnp_common = sum(p$k %in% snp_key))
  drop_key <- snp_key
  EXCL_BED <- Sys.getenv("GT_REDI_EXCLUDE_BED")
  if (nzchar(EXCL_BED)) {
    ex <- fread(EXCL_BED, header = FALSE, select = c(1, 3), col.names = c("chr", "pos"))[chr %in% CHRS]
    ex_key <- unique(key(ex$chr, ex$pos)); rm(ex)
    drop_log$excluded <- sum(p$k %in% ex_key)
    drop_key <- c(drop_key, ex_key); rm(ex_key)
  }
  p <- p[!k %in% drop_key]
  p[, redi := k %in% redi_key]

  ed <- function(a, n) fifelse(n > 0, round(a / pmax(n, 1), 4), 0)
  UV <- list(); multi_n <- c()
  for (v in names(RULES)) {
    r <- RULES[[v]]
    qf <- p$cov_f >= 10 & p$alt_f >= r[["alt"]] & ed(p$alt_f, p$cov_f) >= r[["edit"]]
    qr <- p$cov_r >= 10 & p$alt_r >= r[["alt"]] & ed(p$alt_r, p$cov_r) >= r[["edit"]]
    rev <- qr & (!qf | p$cov_r > p$cov_f)
    w <- which(qf | qr); rv <- rev[w]
    mx <- pmax(fifelse(rv, p$t1_r[w], p$t1_f[w]), fifelse(rv, p$t2_r[w], p$t2_f[w]))
    th <- fifelse(rv, p$t1_r[w] + p$t2_r[w], p$t1_f[w] + p$t2_f[w])
    cv <- fifelse(rv, p$cov_r[w], p$cov_f[w])
    multi <- mx >= 3 & mx >= 0.05 * (cv + th)
    multi_n[[v]] <- sum(multi)
    UV[[v]] <- data.table(k = p$k[w], redi = p$redi[w], orient = fifelse(rv, "reverse", "forward"))[!multi]
    msg(v, ": universe ", nrow(UV[[v]]), " (multiallelic removed ", sum(multi), ")")
  }
  drop_log$multiallelic <- multi_n
  fwrite(rbindlist(lapply(names(drop_log), function(n) {
    x <- unlist(drop_log[[n]]); data.table(removal = n, set = if (is.null(names(x))) "candidates" else names(x), n = x) })),
    file.path(UDIR, "removed_sites.tsv"), sep = "\t")

  gn <- fread(file.path(V9, "truth", "gt_v9_truth_sites.tsv.gz"), select = c("site", "neg_source"))
  gn_key <- unique(site_key(gn[neg_source %in% c("pgp1_ag", "pgp1_hom"), site])); rm(gn)
  msg("PGP1 PASS germline negatives: ", length(gn_key), " (in REDIportal ", sum(gn_key %in% redi_key), ")")

  POS <- lapply(UV, function(u) u[redi == TRUE, k])
  v8p <- fread(file.path(PM, "gt_eval_v8_60/rediportal_pwcov/scoring_mismatch/positives.tsv"))
  pt <- data.table(variant = names(RULES), n_universe = sapply(UV, nrow), n_pos = lengths(POS),
                   n_pos_forward = sapply(UV, function(u) sum(u$redi & u$orient == "forward")),
                   n_pos_reverse = sapply(UV, function(u) sum(u$redi & u$orient == "reverse")),
                   n_neg_universe = sapply(UV, function(u) sum(!u$redi)),
                   n_neg_germline = length(gn_key),
                   n_neg_germline_outside_U = sapply(UV, function(u) sum(!gn_key %in% u$k)))
  pt <- merge(pt, v8p[, .(variant, v8_n_universe = n_universe, v8_n_pos = n_pos)], by = "variant", sort = FALSE)
  fwrite(pt, file.path(UDIR, "positives.tsv"), sep = "\t")
  print(pt)
  for (v in names(POS)) fwrite(data.table(site = key_site(POS[[v]])), file.path(UDIR, paste0("positive_sites_", v, ".tsv.gz")), sep = "\t")
  fwrite(rbindlist(lapply(names(UV), function(v) UV[[v]][redi == TRUE, .(variant = v, site = key_site(k), orient)])),
         file.path(UDIR, "orientation.tsv.gz"), sep = "\t")
  rk <- sort(unique(unlist(lapply(UV, function(u) u[redi == TRUE & orient == "reverse", k]))))
  fwrite(data.table(variant = "any", site = key_site(rk), chr = CHRS[rk %/% 1e9], pos = as.integer(rk %% 1e9)),
         file.path(UDIR, "rev_positive_sites.tsv.gz"), sep = "\t")
  msg("reverse-orientation positive positions (any rule): ", length(rk))
  saveRDS(list(UV = lapply(UV, function(u) u[, .(k, redi, orient)]), POS = POS, gn_key = gn_key, redi_key = redi_key),
          file.path(UDIR, "universe.rds"))
  quit(save = "no")
}

# ---- score -------------------------------------------------------------------------------------------------
OUT <- Sys.getenv("GT_REDI_OUT"); stopifnot(nzchar(OUT)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c("positives.tsv", "removed_sites.tsv", "orientation.tsv.gz", paste0("positive_sites_", names(RULES), ".tsv.gz")))
  file.copy(file.path(UDIR, f), file.path(OUT, f), overwrite = TRUE)
RES <- c(arm = "results_gt_v8", sample = "results_gt_persample_v8", celltype = "results_gt_celltype_v8")
if (nzchar(Sys.getenv("GT_RES_DIRS"))) {
  kv <- strsplit(strsplit(Sys.getenv("GT_RES_DIRS"), ",")[[1]], "=")
  RES[vapply(kv, `[`, "", 1)] <- vapply(kv, `[`, "", 2)
}
CALL_RE <- Sys.getenv("GT_CALL_RE", "_step5_nomm\\.tsv\\.gz$")
MIN_COV <- as.integer(Sys.getenv("GT_UNIT_MIN_COV", "10"))
UNIT_RE <- c(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
TIER_OUT <- c(arm = "dataset", sample = "sample", celltype = "celltype")
UNREPORTED_BAM <- "sample"
UNREPORTED_CALLER_BAM <- list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))

U <- readRDS(file.path(UDIR, "universe.rds"))
POS <- U$POS
files <- rbindlist(lapply(names(RES), function(tr) {
  root <- file.path(PM, RES[[tr]])
  units <- grep(UNIT_RE[[tr]], list.dirs(root, recursive = FALSE, full.names = FALSE), value = TRUE)
  if (nzchar(Sys.getenv("GT_UNIT_KEEP")) && tr == "celltype") {
    keep <- readLines(Sys.getenv("GT_UNIT_KEEP")); stopifnot(all(keep %in% units)); units <- intersect(units, keep)
  }
  rbindlist(lapply(units, function(u) {
    base <- file.path(root, u, "04_filtered")
    f <- list.files(base, CALL_RE, recursive = TRUE, full.names = TRUE)
    p <- tstrsplit(sub(paste0("^", base, "/"), "", f), "/", fixed = TRUE)
    data.table(tier = tr, unit = u, caller = p[[1]], bam_type = p[[2]], strand_mode = p[[3]],
               snp_subtract = sub("^snp_", "", p[[4]]), file = f)
  }))
}))
files <- files[!(bam_type %in% UNREPORTED_BAM |
                 mapply(function(c, b) b %in% UNREPORTED_CALLER_BAM[[c]], caller, bam_type) &
                 caller %in% names(UNREPORTED_CALLER_BAM))]
files[, config := paste(caller, bam_type, strand_mode, snp_subtract, sep = "__")]
msg("reported call files: ", nrow(files))

lk <- data.table(k = unique(c(U$redi_key, U$gn_key, unlist(lapply(U$UV, `[[`, "k")))))
lk[, `:=`(redi_all = k %in% U$redi_key, germ = k %in% U$gn_key)]
for (v in names(RULES)) {
  set(lk, j = paste0("scored_", v), value = lk$k %in% U$UV[[v]]$k | lk$germ)
  set(lk, j = paste0("pos_", v), value = lk$k %in% POS[[v]])
}
setkey(lk, k)
fwd_any <- unique(unlist(lapply(U$UV, function(u) u[redi == TRUE & orient == "forward", k])))
rev_any <- unique(unlist(lapply(U$UV, function(u) u[redi == TRUE & orient == "reverse", k])))
ORI <- lapply(U$UV, function(u) u[redi == TRUE, .(k, orient)])
rm(U)

res <- list()
for (uu in unique(files[, paste(tier, unit)])) {
  tr <- sub(" .*", "", uu); u <- sub("^[^ ]+ ", "", uu)
  pu <- if (u == "NTC_ctrl") "merged_NTC_ctrl" else u
  pdir <- file.path(RD, "pileup", pu)
  stopifnot(length(list.files(pdir, "^chr.*\\.tsv\\.gz$")) == length(CHRS))
  cov_f <- unlist(lapply(CHRS, function(c) {
    p <- fread(file.path(pdir, paste0(c, ".tsv.gz")), select = c("chr", "pos", "coverage"), showProgress = FALSE)[coverage >= MIN_COV]
    kk <- key(p$chr, p$pos); kk[kk %in% fwd_any]
  }))
  rf <- file.path(REVCOV, "ko_truth", paste0(pu, ".tsv.gz"))
  stopifnot(file.exists(rf))
  rp <- fread(rf, select = 1:4, col.names = c("chr", "pos", "ref", "coverage"), showProgress = FALSE)[coverage >= MIN_COV]
  cov_r <- key(rp$chr, rp$pos); cov_r <- cov_r[cov_r %in% rev_any]
  pos_cov <- lapply(ORI, function(o) o[(orient == "forward" & k %in% cov_f) | (orient == "reverse" & k %in% cov_r), k])
  nps_u <- lengths(pos_cov)
  fu <- files[tier == tr & unit == u]
  for (i in seq_len(nrow(fu))) {
    ck <- if (file.size(fu$file[i]) == 0) numeric() else {
      x <- fread(fu$file[i], select = c("site", "varType"), showProgress = FALSE)[varType == "AG"]
      unique(site_key(x$site))
    }
    hit <- lk[J(ck), nomatch = NULL]
    rows <- lapply(names(RULES), function(v) {
      n_sc <- sum(hit[[paste0("scored_", v)]]); tp <- sum(hit[[paste0("pos_", v)]])
      tp_cov <- sum(pos_cov[[v]] %in% ck)
      prec <- if (n_sc) tp / n_sc else NA_real_
      rec <- if (nps_u[[v]]) tp_cov / nps_u[[v]] else NA_real_
      data.table(variant = v, n_calls_in_U = n_sc, n_fp_germline = sum(hit$germ & !hit[[paste0("pos_", v)]]),
                 tp = tp, tp_cov = tp_cov, n_pos = length(POS[[v]]),
                 n_pos_unit = nps_u[[v]], precision = prec, recall = rec,
                 f1 = if (is.na(rec)) NA_real_ else if (is.na(prec) || prec + rec == 0) 0 else 2 * prec * rec / (prec + rec))
    })
    res[[length(res) + 1]] <- cbind(fu[i, .(tier, unit, config, caller, bam_type, strand_mode, snp_subtract)],
                                    n_calls = length(ck), n_known_all = sum(hit$redi_all), rbindlist(rows))
  }
  msg(uu, ": ", nrow(fu), " configurations; covered positives ", paste(names(nps_u), nps_u, collapse = ", "))
}
us <- rbindlist(res)
setcolorder(us, "variant")
fwrite(us, file.path(OUT, "unit_scores.tsv"), sep = "\t")
for (v in names(RULES)) {
  s <- us[variant == v, .(n_discovered = mean(n_calls), n_calls_in_U = mean(n_calls_in_U), n_true = mean(tp),
                          precision = mean(precision, na.rm = TRUE), recall_vs_union = mean(recall, na.rm = TRUE),
                          f1 = mean(f1, na.rm = TRUE), n_true_total = mean(n_pos), n_pos_unit = mean(n_pos_unit),
                          n_units = .N),
          by = .(tier, caller, bam_type, strand_mode, snp_subtract)]
  s[, tier := TIER_OUT[tier]]
  fwrite(s, file.path(OUT, paste0("gt_tier_scoring_", v, ".tsv")), sep = "\t")
  msg(v, ": ", nrow(s), " rows (", uniqueN(s$tier), " tiers)")
}
