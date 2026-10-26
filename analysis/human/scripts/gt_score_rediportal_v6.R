#!/usr/bin/env Rscript
# =============================================================================
# gt_score_rediportal_v6.R   (HANDOFF 2026-09-29 UPDATE-25)
#
# Every reported configuration scored against REDIportal v2 instead of the v6 KO-depletion
# truth, to test whether the ranking depends on how the v6 truth is defined (NTC editing
# thresholds + KO loss). Sites are processed exactly as for the v6 truth
# (run_gt_truth_rebuild_v6.sh), minus the editing-evidence and KO gates:
#   - coverage from the merged NTC BAM (02_gt_bams_umi_v6/merged_NTC_ctrl.bam), samtools-1.20
#     mpileup -B -q30 -Q30 -d 1e6 | parse_pileup_stranded.pl forward (rediportal_pileup_v6.sh;
#     identical to the v6 scan), stranded by the reference base, i.e. REDIportal's strand
#   - dbSNP156 common excluded (gt_eval_v6/snp/dbSNP156_common.bed, as the v6 candidates)
#   - one global positive set, used as the recall denominator in every tier (as v6)
#   - precision only on calls inside the v6 candidate universe U (candidates_noSNP: merged NTC
#     cov >= 5, alt >= 1, edit >= 0.01, not dbSNP common), as v6 scores only labelled sites
#
#   positives (REDIportal, not dbSNP common, merged NTC; nine-truth suite, 2026-10-01):
#       cov10_alt1  cov >= 10, alt >= 1 (edited in the pooled controls)
#       alt3        cov >= 10, alt >= 3
#       alt3_e05    cov >= 10, alt >= 3, edit >= 0.05
#   The v6 candidate floor (cov >= 5, alt >= 1, edit >= 0.01) is computed only as a parity
#   check against candidates_annot (positives.tsv), not scored.
#   call set  = A>G rows of the final call set (<unit>_..._pre100.tsv.gz), chr:pos match, as gt_score_v6.R
#   recall    = positives called / positives
#   precision = calls in U at a REDIportal position / calls in U (REDIportal has no negatives:
#               non-REDIportal candidates count against a method)
#   F1        = 2PR / (P + R); tier value = mean over units (gt_score_v6.R convention)
#
# Inputs : gt_eval_v6/rediportal/{bed/chr*.bed, pileup/merged_NTC_ctrl/chr*.tsv.gz,
#          rediportal_dbsnp_common_sites.txt}; gt_eval_v6/candidates/candidates_annot.tsv
# Overrides (other window states): GT_RES_TEMPLATE (results tree per tier, relative to human/pm,
#   %s = arm/sample/celltype), GT_CALL_RE (call-file regex), GT_REDI_OUT (output dir)
# Optional site removal (removed positions are unscored: out of the positives and of U):
#   GT_REDI_EXCLUDE_BED  BED of positions to drop (e.g. pgp1_snps/PGP1_hu43860C_hg38_snvs.bed; end = pos)
#   GT_REDI_DROP_MULTI   default 1 (user, 2026-10-12; 0 keeps them, as all scorings before that date):
#                        drop multiallelic sites in the merged NTC (gt_multiallelic_counts_v6.R rule:
#                        third allele on the transcribed strand >= 3 reads and >= 5% of depth);
#                        REDIportal positives from rediportal/pileup_full, U from pileups/bq30
# Outputs: gt_eval_v6/rediportal/scoring/
#   positives.tsv                 positives per variant (+ check of floor against candidates_annot)
#   positive_sites_<variant>.tsv.gz  chr:pos of each variant's positives
#   unit_scores.tsv               variant x tier x unit x configuration
#   gt_tier_scoring_<variant>.tsv reported configurations, gt_tier_scoring.tsv format
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("SCRNA_DIR", ".")
PM <- file.path(B, "human/pm")
V6 <- file.path(PM, "gt_eval_v6")
RD <- file.path(V6, "rediportal")
OUT <- Sys.getenv("GT_REDI_OUT", file.path(RD, "scoring"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RES <- c(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6")
if (nzchar(Sys.getenv("GT_RES_TEMPLATE"))) RES[] <- sprintf(Sys.getenv("GT_RES_TEMPLATE"), names(RES))
CALL_RE <- Sys.getenv("GT_CALL_RE", "_pre100\\.tsv\\.gz$")
UNIT_RE <- c(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
TIER_OUT <- c(arm = "dataset", sample = "sample", celltype = "celltype")
UNREPORTED_BAM <- "sample"
UNREPORTED_CALLER_BAM <- list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))
CHRS <- c(paste0("chr", 1:22), "chrX", "chrY")
key <- function(chr, pos) match(chr, CHRS) * 1e9 + as.numeric(pos)
site_key <- function(s) {
  if (!length(s)) return(numeric())
  p <- tstrsplit(s, ":", fixed = TRUE); key(p[[1]], p[[2]])
}
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

# ---- REDIportal positions, dbSNP common removed ----------------------------------
redi <- rbindlist(lapply(CHRS, function(c)
  fread(file.path(RD, "bed", paste0(c, ".bed")), header = FALSE, select = c(1, 3), col.names = c("chr", "pos"))))
redi_key <- unique(key(redi$chr, redi$pos)); rm(redi)
snp_key <- site_key(fread(file.path(RD, "rediportal_dbsnp_common_sites.txt"), header = FALSE)$V1)
redi_nosnp <- setdiff(redi_key, snp_key)
msg("REDIportal positions: ", length(redi_key), "; dbSNP common: ", length(snp_key),
    "; kept: ", length(redi_nosnp))

# ---- v6 candidate universe U -------------------------------------------------------
ca <- fread(file.path(V6, "candidates", "candidates_annot.tsv"))
U <- site_key(ca[in_dbsnp_common == 0, site])
U_redi_v6 <- site_key(ca[in_dbsnp_common == 0 & in_rediportal == 1, site])
msg("v6 candidate universe (no dbSNP common): ", length(U), "; REDIportal among them: ", length(U_redi_v6))

# ---- optional removal: excluded positions, multiallelic sites ----------------------------
drop_key <- numeric(); drop_log <- list()
EXCL_BED <- Sys.getenv("GT_REDI_EXCLUDE_BED")
if (nzchar(EXCL_BED)) {
  ex <- fread(EXCL_BED, header = FALSE, select = c(1, 3), col.names = c("chr", "pos"))[chr %in% CHRS]
  ex_key <- unique(key(ex$chr, ex$pos)); rm(ex)
  drop_log$excluded <- c(redi_nosnp = sum(redi_nosnp %in% ex_key), U = sum(U %in% ex_key))
  msg("excluded positions (", basename(EXCL_BED), "): ", length(ex_key), "; in REDIportal (no SNP): ",
      drop_log$excluded[1], "; in U: ", drop_log$excluded[2])
  drop_key <- c(drop_key, ex_key); rm(ex_key)
}
if (Sys.getenv("GT_REDI_DROP_MULTI", "1") == "1") {
  is_multi <- function(p) {
    isA <- toupper(p$ref) == "A"
    t1 <- fifelse(isA, p$C, p$a); t2 <- fifelse(isA, p$T, p$g); mx <- pmax(t1, t2)
    mx >= 3 & mx >= 0.05 * (p$coverage + t1 + t2)
  }
  pf <- rbindlist(lapply(CHRS, function(c)
    fread(file.path(RD, "pileup_full", "merged_NTC_ctrl", paste0(c, ".tsv.gz")), showProgress = FALSE)[
      coverage >= 10 & variant_cov >= 1]))
  pf[, k := key(chr, pos)]
  pf <- pf[k %in% redi_nosnp]; m_redi <- pf$k[is_multi(pf)]; rm(pf)
  bq <- fread(file.path(V6, "pileups", "bq30", "pileup_merged_NTC_ctrl.tsv"),
              select = c("chr", "pos", "ref", "coverage", "C", "T", "a", "g"), showProgress = FALSE)
  bq[, k := key(chr, pos)]
  bq <- unique(bq[k %in% U], by = "k")
  stopifnot(nrow(bq) == length(U))
  m_U <- bq$k[is_multi(bq)]; rm(bq)
  drop_log$multiallelic <- c(redi_nosnp = length(m_redi), U = length(m_U))
  msg("multiallelic: REDIportal (cov >= 10, alt >= 1) ", length(m_redi), "; U ", length(m_U))
  drop_key <- c(drop_key, m_redi, m_U)
}
if (length(drop_key)) {
  drop_key <- unique(drop_key)
  redi_nosnp <- setdiff(redi_nosnp, drop_key); U <- setdiff(U, drop_key); U_redi_v6 <- setdiff(U_redi_v6, drop_key)
  fwrite(rbindlist(lapply(names(drop_log), function(n) data.table(removal = n, set = names(drop_log[[n]]),
                                                                    n = drop_log[[n]]))),
         file.path(OUT, "removed_sites.tsv"), sep = "\t")
  msg("after removal: REDIportal (no SNP) ", length(redi_nosnp), "; U ", length(U))
}

# ---- positives from the merged NTC pileup --------------------------------------------
pu <- rbindlist(lapply(CHRS, function(c)
  fread(file.path(RD, "pileup", "merged_NTC_ctrl", paste0(c, ".tsv.gz")),
        select = c("chr", "pos", "coverage", "variant_cov"), showProgress = FALSE)))
pu[, k := key(chr, pos)]
pu <- pu[k %in% redi_nosnp]
pu[, edit := round(variant_cov / coverage, 4)]
floor_k <- pu[coverage >= 5 & variant_cov >= 1 & edit >= 0.01, k]
POS <- list(cov10_alt1 = pu[coverage >= 10 & variant_cov >= 1, k],
            alt3 = pu[coverage >= 10 & variant_cov >= 3, k],
            alt3_e05 = pu[coverage >= 10 & variant_cov >= 3 & edit >= 0.05, k])
rm(pu)
chk <- c(floor_in_U_redi = sum(floor_k %in% U_redi_v6), U_redi_in_floor = sum(U_redi_v6 %in% floor_k))
fwrite(data.table(variant = names(POS), n_pos = lengths(POS),
                  n_in_v6_candidates = sapply(POS, function(p) sum(p %in% U)))[
         , `:=`(v6_candidates_redi_nosnp = length(U_redi_v6), n_floor = length(floor_k),
                floor_check_both_ways = paste(chk, collapse = "/"))],
       file.path(OUT, "positives.tsv"), sep = "\t")
for (v in names(POS))
  fwrite(data.table(site = sprintf("%s:%d", CHRS[POS[[v]] %/% 1e9], as.integer(POS[[v]] %% 1e9))),
         file.path(OUT, paste0("positive_sites_", v, ".tsv.gz")), sep = "\t")
msg("positives: ", paste(names(POS), lengths(POS), collapse = ", "),
    "; floor vs v6 REDIportal candidates: ", paste(names(chk), chk, collapse = ", "))
if (chk[1] != length(floor_k) || chk[2] != length(U_redi_v6))
  warning("floor positives differ from the v6 REDIportal candidates; see positives.tsv")
rm(floor_k)

# ---- call files ----------------------------------------------------------------------
files <- rbindlist(lapply(names(RES), function(tr) {
  root <- file.path(PM, RES[[tr]])
  units <- grep(UNIT_RE[[tr]], list.dirs(root, recursive = FALSE, full.names = FALSE), value = TRUE)
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

# ---- lookup table: one keyed row per REDIportal or candidate-universe position ------------
lk <- data.table(k = unique(c(redi_key, U)))
lk[, `:=`(in_U = k %in% U, redi_all = k %in% redi_key, redi_nosnp = k %in% redi_nosnp)]
for (v in names(POS)) set(lk, j = paste0("pos_", v), value = lk$k %in% POS[[v]])
setkey(lk, k)
rm(redi_key, redi_nosnp, U, U_redi_v6)
msg("lookup table: ", nrow(lk), " positions")

# ---- score ---------------------------------------------------------------------------
# Global positives (merged NTC coverage, one denominator for every tier, as v6).
res <- list()
for (uu in unique(files[, paste(tier, unit)])) {
  tr <- sub(" .*", "", uu); u <- sub("^[^ ]+ ", "", uu)
  fu <- files[tier == tr & unit == u]
  for (i in seq_len(nrow(fu))) {
    ck <- if (file.size(fu$file[i]) == 0) numeric() else {
      x <- fread(fu$file[i], select = c("site", "varType"), showProgress = FALSE)[varType == "AG"]
      unique(site_key(x$site))
    }
    hit <- lk[J(ck), nomatch = NULL]
    n_inU <- sum(hit$in_U)
    n_known_U <- sum(hit$in_U & hit$redi_nosnp)
    prec <- if (n_inU) n_known_U / n_inU else NA_real_
    tps <- sapply(names(POS), function(v) sum(hit[[paste0("pos_", v)]]))
    nps <- lengths(POS)
    rec <- ifelse(nps > 0, tps / nps, NA_real_)
    res[[length(res) + 1]] <- data.table(variant = names(tps), fu[i, .(tier, unit, config, caller, bam_type, strand_mode, snp_subtract)],
      n_calls = length(ck), n_calls_in_U = n_inU, n_known_in_U = n_known_U, n_known_all = sum(hit$redi_all),
      n_pos = nps, tp = tps, precision = prec, recall = rec,
      f1 = ifelse(!is.na(prec) & !is.na(rec) & prec + rec > 0, 2 * prec * rec / (prec + rec), NA_real_))
  }
  msg(uu, ": ", nrow(fu), " configurations")
}
us <- rbindlist(res)
fwrite(us, file.path(OUT, "unit_scores.tsv"), sep = "\t")

for (v in names(POS)) {
  s <- us[variant == v, .(n_discovered = mean(n_calls), n_calls_in_U = mean(n_calls_in_U), n_true = mean(tp),
                          precision = mean(precision, na.rm = TRUE), recall_vs_union = mean(recall, na.rm = TRUE),
                          f1 = mean(f1, na.rm = TRUE), n_true_total = mean(n_pos), n_units = .N),
          by = .(tier, caller, bam_type, strand_mode, snp_subtract)]
  s[, tier := TIER_OUT[tier]]
  fwrite(s, file.path(OUT, paste0("gt_tier_scoring_", v, ".tsv")), sep = "\t")
  msg(v, ": ", nrow(s), " rows (", uniqueN(s$tier), " tiers)")
}
