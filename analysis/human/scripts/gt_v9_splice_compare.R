#!/usr/bin/env Rscript
# gt_v9_splice_compare.R splice|compare   (2026-11-07)
#   splice   GT_V9_DIR/_parts/{pub,ch}/<state>/<tier>/ (58 published + 2 Cellsnp-lite hybrid configurations, as
#            gt_v8_splice60.R) -> GT_V9_DIR/scoring_pwcov/<state>/<tier>/
#   splice76 scoring_pwcov/<state>/<tier>/ (60) + umimd/<tier>/ (16 UMI-aware MarkDuplicates configurations,
#            umimd_score.sh) -> GT_V9_DIR/scoring_pwcov76/<state>/<tier>/   (2026-11-13)
#   compare  v9 ranks vs gt_eval_v8_60 (same state) -> GT_V9_DIR/compare_v8_<state>.tsv + printed summary
suppressPackageStartupMessages(library(data.table))
PM = file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
V9 = Sys.getenv("GT_V9_DIR", file.path(PM, "gt_eval_v9"))
STATE = Sys.getenv("GT_STATE", "mismatch")
CUR = file.path(PM, "gt_eval_v8_60/scoring_pwcov", STATE)
O = file.path(V9, "scoring_pwcov", STATE)
TIERS = c("arm", "sample", "celltype")
FILES = c("gt_v7_scores.tsv", "gt_v7_scores_reported.tsv", "gt_v7_unit_scores.tsv", "gt_v7_calls_NTC.tsv.gz")
part = commandArgs(trailingOnly = TRUE)[1]

if (part == "splice") {
  for (tr in TIERS) for (f in FILES) {
    x = fread(file.path(V9, "_parts/pub", STATE, tr, f)); y = fread(file.path(V9, "_parts/ch", STATE, tr, f))
    stopifnot(identical(names(x), names(y)), !any(unique(y$config) %in% x$config), uniqueN(y$config) == 2L)
    if (f == "gt_v7_scores_reported.tsv") stopifnot(uniqueN(x$config) == 58L)
    dir.create(file.path(O, tr), recursive = TRUE, showWarnings = FALSE)
    fwrite(rbind(x, y), file.path(O, tr, f), sep = "\t")
  }
  cat("spliced", STATE, "\n")
}

if (part == "splice76") {
  O76 = file.path(V9, "scoring_pwcov76", STATE)
  for (tr in TIERS) for (f in FILES) {
    x = fread(file.path(O, tr, f)); y = fread(file.path(V9, "umimd", tr, f))
    stopifnot(identical(names(x), names(y)), !any(unique(y$config) %in% x$config), uniqueN(y$config) == 16L,
              all(y$bam_type == "sampleUmiMd"))
    if (f == "gt_v7_scores_reported.tsv") stopifnot(uniqueN(x$config) == 60L)
    dir.create(file.path(O76, tr), recursive = TRUE, showWarnings = FALSE)
    fwrite(rbind(x, y), file.path(O76, tr, f), sep = "\t")
  }
  cat("spliced76", STATE, "\n")
}

if (part == "compare") {
  K = c("gate", "t", "config", "caller", "bam_type", "strand_mode", "snp_subtract")
  V = c("overall_rank", "mean_f1", paste0("f1_", TIERS), paste0("precision_", TIERS), paste0("recall_", TIERS))
  a = fread(file.path(CUR, "config_rank.tsv"))[, c(K, V), with = FALSE]
  b = fread(file.path(O, "config_rank.tsv"))[, c(K, V), with = FALSE]
  m = merge(a, b, by = K, suffixes = c("_v8", "_v9"))
  stopifnot(nrow(m) == nrow(a), nrow(a) == nrow(b), m[, .N, by = .(gate, t)][, all(N == 60L)])
  m[, d_rank := overall_rank_v9 - overall_rank_v8]
  fwrite(m, file.path(V9, paste0("compare_v8_", STATE, ".tsv")), sep = "\t")
  cat("\n== rank agreement v8 vs v9 ==\n")
  print(m[, .(spearman = round(cor(overall_rank_v8, overall_rank_v9, method = "spearman"), 3),
              max_abs_shift = max(abs(d_rank)), top1_v8 = config[overall_rank_v8 == 1][1],
              top1_v9 = config[overall_rank_v9 == 1][1]), by = .(gate, t)])
  cat("\n== t50: by stranding x SNP removal (mean rank, mean F1 change, mean precision per tier under v9) ==\n")
  print(m[t == "t50", .(n = .N, rank_v8 = round(mean(overall_rank_v8), 1), rank_v9 = round(mean(overall_rank_v9), 1),
                        d_f1 = signif(mean(mean_f1_v9 - mean_f1_v8), 3),
                        prec_arm = round(mean(precision_arm_v9), 3), prec_sample = round(mean(precision_sample_v9), 3),
                        prec_ct = round(mean(precision_celltype_v9), 3)), by = .(gate, strand_mode, snp_subtract)][order(gate, rank_v9)])
  for (g in c("depth", "expected")) {
    cat("\n== t50", g, ": top 12 under v9 ==\n")
    print(m[t == "t50" & gate == g][order(overall_rank_v9)][1:12, .(config, rank_v8 = overall_rank_v8, rank_v9 = overall_rank_v9,
          f1_v9 = round(mean_f1_v9, 3), f1_arm = round(f1_arm_v9, 3), f1_sample = round(f1_sample_v9, 3),
          f1_ct = round(f1_celltype_v9, 3), prec_arm = round(precision_arm_v9, 3))])
  }
  cr = function(s) m[t == "t50", .(rank = min(get(paste0("overall_rank_", s)))), by = .(gate, caller)]
  cat("\n== t50: caller rank (best configuration) ==\n")
  print(merge(cr("v8"), cr("v9"), by = c("gate", "caller"), suffixes = c("_v8", "_v9"))[order(gate, rank_v9)])
}
