#!/usr/bin/env Rscript
# =============================================================================
# gt_score_rediportal_v8.R   (2026-10-13)
#
# REDIportal v2 as truth, v8 rules (user 2026-10-13). Differs from gt_score_rediportal_v6.R in two ways:
#   1. per-truth universe U_v: merged NTC A-to-G candidate positions (gt_eval_v6/pileups/bq30, the v6 scan:
#      stranded by the reference base, BQ/MQ >= 30) that pass the truth's own site rule plus an editing
#      level floor of 1%; dbSNP common, GT_REDI_EXCLUDE_BED (PGP1) and multiallelic positions removed.
#        cov10_alt1  cov >= 10, alt >= 1, edit >= 0.01
#        alt3        cov >= 10, alt >= 3, edit >= 0.01
#        alt3_e05    cov >= 10, alt >= 3, edit >= 0.05
#      positives_v = U_v at a REDIportal position; negatives_v = the rest of U_v
#   2. per-unit recall denominator: positives_v with >= GT_UNIT_MIN_COV (10) stranded reads in the unit
#      (gt_eval_v6/rediportal/pileup/<unit>/, rediportal_pileup_v6.sh; the arm unit is merged_NTC_ctrl)
#   call set  = A>G rows of the call files, chr:pos match (as v6)
#   precision = calls at positives_v / calls in U_v
#   recall    = calls at covered positives_v / covered positives_v (NA if none covered)
#   F1        = 2PR / (P + R); 0 when the unit has covered positives but no hit (including no calls in U_v);
#               NA only without covered positives. Tier F1 = mean of the unit F1 (f1 column; units as gt_score_v7.R)
# Positives are checked against the REDIportal pileup of the merged NTC (positives.tsv).
#
# Overrides: GT_RES_TEMPLATE, GT_CALL_RE, GT_REDI_OUT, GT_REDI_EXCLUDE_BED, GT_REDI_DROP_MULTI (as v6),
#            GT_UNIT_MIN_COV (default 10), GT_UNIT_KEEP (file, one cell-type unit per line: celltype tier restricted to them)
# Outputs (GT_REDI_OUT): positives.tsv, positive_sites_<variant>.tsv.gz, removed_sites.tsv, unit_scores.tsv,
#   gt_tier_scoring_<variant>.tsv (gt_tier_scoring.tsv format, as v6)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("SCRNA_DIR", ".")
PM <- file.path(B, "human/pm")
V6 <- file.path(PM, "gt_eval_v6")
RD <- file.path(V6, "rediportal")
OUT <- Sys.getenv("GT_REDI_OUT", file.path(PM, "gt_eval_v8", "rediportal", "scoring_mismatch"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RES <- c(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6")
if (nzchar(Sys.getenv("GT_RES_TEMPLATE"))) RES[] <- sprintf(Sys.getenv("GT_RES_TEMPLATE"), names(RES))
# GT_RES_DIRS="arm=<dir>,sample=<dir>,celltype=<dir>" (relative to human/pm) for result trees not named by tier
if (nzchar(Sys.getenv("GT_RES_DIRS"))) {
  kv <- strsplit(strsplit(Sys.getenv("GT_RES_DIRS"), ",")[[1]], "=")
  RES[vapply(kv, `[`, "", 1)] <- vapply(kv, `[`, "", 2)
}
CALL_RE <- Sys.getenv("GT_CALL_RE", "_pre100\\.tsv\\.gz$")
MIN_COV <- as.integer(Sys.getenv("GT_UNIT_MIN_COV", "10"))
UNIT_RE <- c(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
TIER_OUT <- c(arm = "dataset", sample = "sample", celltype = "celltype")
UNREPORTED_BAM <- "sample"
UNREPORTED_CALLER_BAM <- list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))
CHRS <- c(paste0("chr", 1:22), "chrX", "chrY")
RULES <- list(cov10_alt1 = c(alt = 1, edit = 0.01), alt3 = c(alt = 3, edit = 0.01), alt3_e05 = c(alt = 3, edit = 0.05))
key <- function(chr, pos) match(chr, CHRS) * 1e9 + as.numeric(pos)
site_key <- function(s) {
  if (!length(s)) return(numeric())
  p <- tstrsplit(s, ":", fixed = TRUE); key(p[[1]], p[[2]])
}
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

# ---- REDIportal positions, dbSNP common ---------------------------------------------
redi <- rbindlist(lapply(CHRS, function(c)
  fread(file.path(RD, "bed", paste0(c, ".bed")), header = FALSE, select = c(1, 3), col.names = c("chr", "pos"))))
redi_key <- unique(key(redi$chr, redi$pos)); rm(redi)
snp_key <- unique(c(site_key(fread(file.path(RD, "rediportal_dbsnp_common_sites.txt"), header = FALSE)$V1),
                    site_key(fread(file.path(V6, "candidates", "candidates_annot.tsv"))[in_dbsnp_common == 1, site])))
msg("REDIportal positions: ", length(redi_key), "; dbSNP common (REDIportal + candidates): ", length(snp_key))

# ---- merged NTC candidate pileup (the universe source) --------------------------------
bq <- fread(file.path(V6, "pileups", "bq30", "pileup_merged_NTC_ctrl.tsv"),
            select = c("chr", "pos", "ref", "coverage", "variant_cov", "C", "T", "a", "g"), showProgress = FALSE)
bq[, k := key(chr, pos)]
bq <- unique(bq, by = "k")
bq[, edit := round(variant_cov / coverage, 4)]
msg("merged NTC candidate positions: ", nrow(bq))

# ---- removal: dbSNP common, excluded positions, multiallelic ---------------------------
drop_log <- list(dbsnp_common = c(candidates = sum(bq$k %in% snp_key)))
drop_key <- snp_key
EXCL_BED <- Sys.getenv("GT_REDI_EXCLUDE_BED")
if (nzchar(EXCL_BED)) {
  ex <- fread(EXCL_BED, header = FALSE, select = c(1, 3), col.names = c("chr", "pos"))[chr %in% CHRS]
  ex_key <- unique(key(ex$chr, ex$pos)); rm(ex)
  drop_log$excluded <- c(candidates = sum(bq$k %in% ex_key))
  msg("excluded positions (", basename(EXCL_BED), "): ", length(ex_key), "; among candidates: ", drop_log$excluded)
  drop_key <- c(drop_key, ex_key); rm(ex_key)
}
if (Sys.getenv("GT_REDI_DROP_MULTI", "1") == "1") {
  isA <- toupper(bq$ref) == "A"
  t1 <- fifelse(isA, bq$C, bq$a); t2 <- fifelse(isA, bq$T, bq$g); mx <- pmax(t1, t2)
  m_k <- bq$k[mx >= 3 & mx >= 0.05 * (bq$coverage + t1 + t2)]
  drop_log$multiallelic <- c(candidates = length(m_k))
  msg("multiallelic candidates: ", length(m_k))
  drop_key <- c(drop_key, m_k); rm(isA, t1, t2, mx, m_k)
}
drop_key <- unique(drop_key)
bq <- bq[!k %in% drop_key]
fwrite(rbindlist(lapply(names(drop_log), function(n) data.table(removal = n, set = names(drop_log[[n]]), n = drop_log[[n]]))),
       file.path(OUT, "removed_sites.tsv"), sep = "\t")
bq[, redi := k %in% redi_key]

# ---- per-truth universes and positives --------------------------------------------------
UV <- lapply(RULES, function(r) bq[coverage >= 10 & variant_cov >= r[["alt"]] & edit >= r[["edit"]], .(k, redi)])
POS <- lapply(UV, function(u) u[redi == TRUE, k])

# parity: positives from the REDIportal pileup of the merged NTC under the same rules
pu <- rbindlist(lapply(CHRS, function(c)
  fread(file.path(RD, "pileup", "merged_NTC_ctrl", paste0(c, ".tsv.gz")),
        select = c("chr", "pos", "coverage", "variant_cov"), showProgress = FALSE)[coverage >= 10 & variant_cov >= 1]))
pu[, k := key(chr, pos)]
pu <- pu[!k %in% drop_key][, edit := round(variant_cov / coverage, 4)]
POS_PU <- lapply(RULES, function(r) pu[variant_cov >= r[["alt"]] & edit >= r[["edit"]], k])
rm(pu)
pt <- data.table(variant = names(RULES), n_universe = sapply(UV, nrow), n_pos = lengths(POS),
                 n_neg = sapply(UV, function(u) sum(!u$redi)),
                 n_pos_redi_pileup = lengths(POS_PU),
                 n_pos_both = mapply(function(a, b) sum(a %in% b), POS, POS_PU))
fwrite(pt, file.path(OUT, "positives.tsv"), sep = "\t")
print(pt)
if (any(pt$n_pos_both != pt$n_pos | pt$n_pos_both != pt$n_pos_redi_pileup))
  warning("candidate-pileup positives differ from the REDIportal-pileup positives; see positives.tsv")
rm(POS_PU)
for (v in names(POS))
  fwrite(data.table(site = sprintf("%s:%d", CHRS[POS[[v]] %/% 1e9], as.integer(POS[[v]] %% 1e9))),
         file.path(OUT, paste0("positive_sites_", v, ".tsv.gz")), sep = "\t")

# ---- call files ------------------------------------------------------------------------
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

lk <- data.table(k = unique(c(redi_key, bq$k)))
lk[, redi_all := k %in% redi_key]
for (v in names(RULES)) {
  set(lk, j = paste0("inU_", v), value = lk$k %in% UV[[v]]$k)
  set(lk, j = paste0("pos_", v), value = lk$k %in% POS[[v]])
}
setkey(lk, k)
pos_any <- unique(unlist(POS))
rm(redi_key, bq, UV)

# ---- score -------------------------------------------------------------------------------
res <- list()
for (uu in unique(files[, paste(tier, unit)])) {
  tr <- sub(" .*", "", uu); u <- sub("^[^ ]+ ", "", uu)
  pdir <- file.path(RD, "pileup", if (u == "NTC_ctrl") "merged_NTC_ctrl" else u)
  stopifnot(length(list.files(pdir, "^chr.*\\.tsv\\.gz$")) == length(CHRS))
  cov_k <- unlist(lapply(CHRS, function(c) {
    p <- fread(file.path(pdir, paste0(c, ".tsv.gz")), select = c("chr", "pos", "coverage"), showProgress = FALSE)[coverage >= MIN_COV]
    kk <- key(p$chr, p$pos); kk[kk %in% pos_any]
  }))
  pos_cov <- lapply(POS, function(p) p[p %in% cov_k])
  nps_u <- lengths(pos_cov)
  fu <- files[tier == tr & unit == u]
  for (i in seq_len(nrow(fu))) {
    ck <- if (file.size(fu$file[i]) == 0) numeric() else {
      x <- fread(fu$file[i], select = c("site", "varType"), showProgress = FALSE)[varType == "AG"]
      unique(site_key(x$site))
    }
    hit <- lk[J(ck), nomatch = NULL]
    ck_cov <- ck[ck %in% cov_k]
    rows <- lapply(names(RULES), function(v) {
      n_inU <- sum(hit[[paste0("inU_", v)]]); tp <- sum(hit[[paste0("pos_", v)]])
      tp_cov <- sum(pos_cov[[v]] %in% ck_cov)
      prec <- if (n_inU) tp / n_inU else NA_real_
      rec <- if (nps_u[[v]]) tp_cov / nps_u[[v]] else NA_real_
      data.table(variant = v, n_calls_in_U = n_inU, tp = tp, tp_cov = tp_cov, n_pos = length(POS[[v]]),
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
