#!/usr/bin/env Rscript
# =============================================================================
# gt_score_v6.R
#
# Score every method configuration against the frozen v6 ground truth
# (human/pm/gt_eval_v6/truth/gt_truth_sites.tsv.gz). Discoveries are the NTC
# control's calls only (user, 2026-09-28): the truth is defined by editing in
# NTC plus loss in the KO arms, so a method is scored on the same arm.
#
#   configuration = caller x bam_type x strand_mode x snp_subtract (never pooled)
#   call set      = A>G rows (varType == "AG") of <unit>_..._pre100.tsv.gz,
#                   matched to the truth on chr:pos (strand not used)
#   precision     = TRUE / (TRUE + FP_known_not_depleted + FP_novel);
#                   ambiguous and not-in-truth calls are excluded (counted only)
#   recall        = TRUE called / all ground-truth sites (1,692; the truth is
#                   caller-independent, so the denominator is global)
#   F1            = 2PR / (P + R)
#   reported      = FALSE for the Unprocessed BAM (bam_type "sample"): it is
#                   scored and kept in the full tables but is not reported in
#                   the main results, figures or plots (user, 2026-09-28)
#
# Tiers (one NTC unit per replicate):
#   arm       results_gt_v6/NTC_ctrl                 1 unit; SE by site bootstrap
#   sample    results_gt_persample_v6/T{1..4}_NTC_ctrl  4 units; mean and SE across teratomas
#   celltype  results_gt_celltype_v6/NTC_ctrl__*     n units; mean and SD across cell types
#
# Usage: Rscript gt_score_v6.R <tier> [results_dir] [outdir] [call_file_regex]
#   call_file_regex defaults to _pre100\.tsv\.gz$; the pre-window comparison passes
#   _step4_cov_edit_filtered\.tsv\.gz$ (gt_eval_v6/prewindow/, HANDOFF 2026-09-29 UPDATE-21)
# Outputs (<outdir>, default human/pm/gt_eval_v6/scoring/<tier>/):
#   gt_v6_unit_scores.tsv    one row per unit x configuration
#   gt_v6_scores.tsv         one row per configuration (tier summary; all BAM types)
#   gt_v6_scores_reported.tsv  the same, reported configurations only (no Unprocessed)
#   gt_v6_calls_NTC.tsv.gz   long table of scored calls (unit, configuration, site, class)
# =============================================================================

suppressPackageStartupMessages(library(data.table))
args = commandArgs(trailingOnly = TRUE)
TIER = args[1]
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
DEF = list(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6")
UNIT_RE = list(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
stopifnot(TIER %in% names(DEF))
RES = if (length(args) >= 2) args[2] else file.path(PM, DEF[[TIER]])
OUT = if (length(args) >= 3) args[3] else file.path(PM, "gt_eval_v6", "scoring", TIER)
CALL_RE = if (length(args) >= 4) args[4] else "_pre100\\.tsv\\.gz$"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_BOOT = 500
UNREPORTED_BAM = "sample"   # Unprocessed: recorded, not reported
# A BAM-processing configuration whose defining parameter has no effect on a caller's calls is
# recorded but not reported for that caller (HANDOFF 2026-09-29 UPDATE-19). scAllele has no per-base
# quality filter: recalq20 == recal, and BQSR changes only FILTER, so its recal BAM differs from
# sampleRmDup only by SplitNCigarReads.
UNREPORTED_CALLER_BAM = list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))
unreported = function(caller, bam_type) bam_type %in% UNREPORTED_BAM |
  mapply(function(c, b) b %in% UNREPORTED_CALLER_BAM[[c]], caller, bam_type) &
  caller %in% names(UNREPORTED_CALLER_BAM)
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

# ---- truth ---------------------------------------------------------------------
truth = fread(file.path(PM, "gt_eval_v6", "truth", "gt_truth_sites.tsv.gz"), select = c("site", "class"))
N_TRUE = truth[class == "TRUE", .N]
msg("truth: ", nrow(truth), " sites; ", N_TRUE, " ground-truth (TRUE)")

# ---- call sets -----------------------------------------------------------------
units = grep(UNIT_RE[[TIER]], list.dirs(RES, recursive = FALSE, full.names = FALSE), value = TRUE)
if (!length(units)) stop("no NTC units matching ", UNIT_RE[[TIER]], " in ", RES)
msg("tier ", TIER, ": ", length(units), " NTC unit(s) in ", RES)
files = rbindlist(lapply(units, function(u) {
  f = list.files(file.path(RES, u, "04_filtered"), CALL_RE, recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NULL)
  p = tstrsplit(sub(paste0("^", file.path(RES, u, "04_filtered"), "/"), "", f), "/", fixed = TRUE)
  data.table(unit = u, caller = p[[1]], bam_type = p[[2]], strand_mode = p[[3]],
             snp_subtract = sub("^snp_", "", p[[4]]), file = f)
}))
files[, config := paste(caller, bam_type, strand_mode, snp_subtract, sep = "__")]
msg("call files: ", nrow(files), " (", uniqueN(files$config), " configurations)")
calls = rbindlist(lapply(seq_len(nrow(files)), function(i) {
  if (file.size(files$file[i]) == 0) return(NULL)
  x = fread(files$file[i], select = c("site", "varType"), showProgress = FALSE)
  unique(x[varType == "AG", .(site)])[, `:=`(unit = files$unit[i], config = files$config[i])]
}))
calls = merge(calls, truth, by = "site", all.x = TRUE)
calls[is.na(class), class := "not_in_truth"]
fwrite(calls[, .(unit, config, site, class)], file.path(OUT, "gt_v6_calls_NTC.tsv.gz"), sep = "\t")

# ---- per unit x configuration --------------------------------------------------
score = function(cl) {
  tp = sum(cl == "TRUE"); fk = sum(cl == "FP_known_not_depleted"); fn = sum(cl == "FP_novel")
  p = if (tp + fk + fn > 0) tp / (tp + fk + fn) else NA_real_
  r = tp / N_TRUE
  list(n_calls = length(cl), n_true = tp, n_fp_known = fk, n_fp_novel = fn,
       n_ambiguous = sum(cl == "ambiguous"), n_not_in_truth = sum(cl == "not_in_truth"),
       precision = p, recall = r, f1 = if (!is.na(p) && p + r > 0) 2 * p * r / (p + r) else NA_real_)
}
us = calls[, score(class), by = .(unit, config)]
# configurations with a call file but no A>G calls still score (0 calls)
us = merge(unique(files[, .(unit, config, caller, bam_type, strand_mode, snp_subtract)]), us,
           by = c("unit", "config"), all.x = TRUE)
for (k in c("n_calls", "n_true", "n_fp_known", "n_fp_novel", "n_ambiguous", "n_not_in_truth", "recall"))
  set(us, which(is.na(us[[k]])), k, 0)
setorder(us, config, unit)
fwrite(us, file.path(OUT, "gt_v6_unit_scores.tsv"), sep = "\t")

# ---- tier summary per configuration ---------------------------------------------
se = function(x) { x = x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
sdv = function(x) { x = x[!is.na(x)]; if (length(x) > 1) sd(x) else NA_real_ }
spread = if (TIER == "celltype") sdv else se
summ = us[, .(n_units = .N, n_calls = mean(n_calls), n_true = mean(n_true),
              n_fp = mean(n_fp_known + n_fp_novel), n_ambiguous = mean(n_ambiguous),
              precision = mean(precision, na.rm = TRUE), recall = mean(recall), f1 = mean(f1, na.rm = TRUE),
              precision_err = spread(precision), recall_err = spread(recall), f1_err = spread(f1)),
          by = .(config, caller, bam_type, strand_mode, snp_subtract)]
summ[, err_type := fifelse(TIER == "celltype", "SD across cell types",
                    fifelse(TIER == "sample", "SE across teratomas", "site bootstrap SE"))]

# arm tier: one unit, so the uncertainty is a bootstrap over the truth-labelled
# (TRUE + FP) sites; ambiguous / not-in-truth calls do not enter P, R or F1
if (TIER == "arm") {
  set.seed(1)
  lab = truth[class %in% c("TRUE", "FP_known_not_depleted", "FP_novel")]
  called = calls[class %in% c("TRUE", "FP_known_not_depleted", "FP_novel"), .(config, site)]
  cfgs = summ$config
  hit = sapply(cfgs, function(cf) lab$site %in% called[config == cf, site])   # sites x configs
  is_t = lab$class == "TRUE"
  bs = replicate(N_BOOT, {
    i = sample.int(nrow(lab), replace = TRUE)
    t_i = is_t[i]; h = hit[i, , drop = FALSE]
    tp = colSums(h & t_i); fp = colSums(h & !t_i); nt = sum(t_i)
    p = ifelse(tp + fp > 0, tp / (tp + fp), NA); r = tp / nt
    rbind(p, r, ifelse(!is.na(p) & p + r > 0, 2 * p * r / (p + r), NA))
  })   # 3 x configs x B
  summ[, `:=`(precision_err = apply(bs[1, , ], 1, sd, na.rm = TRUE),
              recall_err = apply(bs[2, , ], 1, sd, na.rm = TRUE),
              f1_err = apply(bs[3, , ], 1, sd, na.rm = TRUE))]
}
summ[, reported := !unreported(caller, bam_type)]
setorder(summ, -f1)
fwrite(summ, file.path(OUT, "gt_v6_scores.tsv"), sep = "\t")
fwrite(summ[reported == TRUE], file.path(OUT, "gt_v6_scores_reported.tsv"), sep = "\t")
msg("wrote ", nrow(summ), " configurations (", summ[reported == TRUE, .N], " reported) -> ", OUT)
print(head(summ[reported == TRUE, .(config, n_units, n_calls = round(n_calls), n_true = round(n_true, 1),
                    precision = round(precision, 3), recall = round(recall, 3), f1 = round(f1, 3))], 10))
