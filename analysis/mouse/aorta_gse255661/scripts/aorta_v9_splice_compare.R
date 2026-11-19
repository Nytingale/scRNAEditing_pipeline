#!/usr/bin/env Rscript
# aorta_v9_splice_compare.R splice|compare   (user 2026-11-07)
#   splice   58 published + 2 Cellsnp-lite hybrid configurations (arm, sample), mismatch state:
#              v9  gt_eval_v9/_parts/{pub,ch}/mismatch/<tier>  -> gt_eval_v9/scoring_pwcov/mismatch/<tier>
#              v8  gt_eval_v8/scoring_pwcov/mismatch + gt_eval_v8/cellsnp_hybrid/scoring_pwcov
#                  -> gt_eval_v9/v8_60/scoring_pwcov/mismatch/<tier>   (v8 rules at 60, for the comparison)
#   compare  v8 (60) vs v9 ranks -> gt_eval_v9/compare_v8_mismatch.tsv + printed summary
#   splice76 these 60 + 16 UMI-aware MarkDuplicates configurations (gt_eval_v9/umimd/<tier>, aorta_umimd_score.sh)
#            -> gt_eval_v9/scoring_pwcov76/mismatch/<tier>   (2026-11-14)
suppressPackageStartupMessages(library(data.table))
A = file.path(Sys.getenv("SCRNA_DIR", "."), "mouse/aorta_gse255661")
V9 = file.path(A, "gt_eval_v9"); V8 = file.path(A, "gt_eval_v8")
TIERS = c("arm", "sample")
FILES = c("gt_v7_scores.tsv", "gt_v7_scores_reported.tsv", "gt_v7_unit_scores.tsv", "gt_v7_calls_NTC.tsv.gz")
part = commandArgs(trailingOnly = TRUE)[1]

splice = function(pub, ch, out) for (tr in TIERS) for (f in FILES) {
  x = fread(file.path(pub, tr, f)); y = fread(file.path(ch, tr, f))
  if (!identical(names(x), names(y))) y = y[, names(x), with = FALSE]
  stopifnot(!any(unique(y$config) %in% x$config), uniqueN(y$config) == 2L)
  if (f == "gt_v7_scores_reported.tsv") stopifnot(uniqueN(x$config) == 58L)
  dir.create(file.path(out, tr), recursive = TRUE, showWarnings = FALSE)
  fwrite(rbind(x, y), file.path(out, tr, f), sep = "\t")
}

if (part == "splice") {
  splice(file.path(V9, "_parts/pub/mismatch"), file.path(V9, "_parts/ch/mismatch"), file.path(V9, "scoring_pwcov/mismatch"))
  splice(file.path(V8, "scoring_pwcov/mismatch"), file.path(V8, "cellsnp_hybrid/scoring_pwcov"),
         file.path(V9, "v8_60/scoring_pwcov/mismatch"))
  cat("spliced\n")
}

if (part == "splice76") {
  O76 = file.path(V9, "scoring_pwcov76/mismatch")
  for (tr in TIERS) for (f in FILES) {
    x = fread(file.path(V9, "scoring_pwcov/mismatch", tr, f)); y = fread(file.path(V9, "umimd", tr, f))
    if (!identical(names(x), names(y))) y = y[, names(x), with = FALSE]
    stopifnot(!any(unique(y$config) %in% x$config), uniqueN(y$config) == 16L, all(y$bam_type == "sampleUmiMd"))
    if (f == "gt_v7_scores_reported.tsv") stopifnot(uniqueN(x$config) == 60L)
    dir.create(file.path(O76, tr), recursive = TRUE, showWarnings = FALSE)
    fwrite(rbind(x, y), file.path(O76, tr, f), sep = "\t")
  }
  cat("spliced 76\n")
}

if (part == "compare") {
  K = c("gate", "t", "config", "caller", "bam_type", "strand_mode", "snp_subtract")
  V = c("overall_rank", "mean_f1", paste0("f1_", TIERS), paste0("precision_", TIERS), paste0("recall_", TIERS))
  a = fread(file.path(V9, "v8_60/scoring_pwcov/mismatch/config_rank.tsv"))[, c(K, V), with = FALSE]
  b = fread(file.path(V9, "scoring_pwcov/mismatch/config_rank.tsv"))[, c(K, V), with = FALSE]
  m = merge(a, b, by = K, suffixes = c("_v8", "_v9"))
  stopifnot(nrow(m) == nrow(a), nrow(a) == nrow(b), m[, .N, by = .(gate, t)][, all(N == 60L)])
  m[, d_rank := overall_rank_v9 - overall_rank_v8]
  fwrite(m, file.path(V9, "compare_v8_mismatch.tsv"), sep = "\t")
  cat("\n== rank agreement v8 vs v9 (60 configurations) ==\n")
  print(m[, .(spearman = round(cor(overall_rank_v8, overall_rank_v9, method = "spearman"), 3),
              max_abs_shift = max(abs(d_rank)), top1_v8 = config[overall_rank_v8 == 1][1],
              top1_v9 = config[overall_rank_v9 == 1][1]), by = .(gate, t)])
  cat("\n== t50: by stranding x SNP removal ==\n")
  print(m[t == "t50", .(n = .N, rank_v8 = round(mean(overall_rank_v8), 1), rank_v9 = round(mean(overall_rank_v9), 1),
                        d_f1 = signif(mean(mean_f1_v9 - mean_f1_v8), 3),
                        prec_arm_v8 = round(mean(precision_arm_v8), 3), prec_arm_v9 = round(mean(precision_arm_v9), 3)),
          by = .(gate, strand_mode, snp_subtract)][order(gate, rank_v9)])
  for (g in c("depth", "expected")) {
    cat("\n== t50", g, ": top 12 under v9 ==\n")
    print(m[t == "t50" & gate == g][order(overall_rank_v9)][1:12, .(config, rank_v8 = overall_rank_v8, rank_v9 = overall_rank_v9,
          f1_v9 = round(mean_f1_v9, 3), prec_arm = round(precision_arm_v9, 3), rec_arm = round(recall_arm_v9, 3))])
  }
  cr = function(s) m[t == "t50", .(rank = min(get(paste0("overall_rank_", s)))), by = .(gate, caller)]
  cat("\n== t50: caller rank (best configuration) ==\n")
  print(merge(cr("v8"), cr("v9"), by = c("gate", "caller"), suffixes = c("_v8", "_v9"))[order(gate, rank_v9)])
}
