#!/usr/bin/env Rscript
# =============================================================================
# gt_rank_v7.R   (2026-10-08)
# Cross-tier ranks on the v7 truth (gt_score_v7.R outputs, reported configurations), per gate x t:
# F1 rank within each tier (ties min); overall rank by mean F1 across arm / sample / celltype (ties: lower rank sum
# across tiers); caller rank from each caller's best configuration
# (as gt_truth_agtc_test_v6.R / TERATOMA_window_state_ranksum_plot.R).
# Usage : Rscript gt_rank_v7.R [scoring_dir]   (default gt_eval_v7/scoring)
# Output: <scoring_dir>/config_rank.tsv, caller_rank.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
args = commandArgs(trailingOnly = TRUE)
SC = if (length(args) >= 1) args[1] else file.path(B, "human/pm/gt_eval_v7/scoring")
TIERS = c("arm", "sample", "celltype")
COMBO = c("caller", "bam_type", "strand_mode", "snp_subtract")
KEY = c("gate", "t")
s = rbindlist(lapply(TIERS, function(tr) fread(file.path(SC, tr, "gt_v7_scores_reported.tsv"))[, tier := tr]))
stopifnot(s[, uniqueN(config), by = c(KEY, "tier")][, uniqueN(V1) == 1])
s[, f1_rank := frank(-fcoalesce(f1, -Inf), ties.method = "min"), by = c(KEY, "tier")]
w = dcast(s, as.formula(paste(paste(c(KEY, "config", COMBO), collapse = " + "), "~ tier")),
          value.var = c("f1", "recall", "precision", "f1_rank"))
rs = s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = c(KEY, "config", COMBO)]
rs[, overall_rank := frank(list(-mean_f1, rank_sum), ties.method = "min"), by = KEY]
rs = merge(rs, w, by = c(KEY, "config", COMBO))
setorderv(rs, c(KEY, "overall_rank"))
fwrite(rs, file.path(SC, "config_rank.tsv"), sep = "\t")
cb = rs[, .SD[1], by = c(KEY, "caller")]
cb[, caller_rank := frank(overall_rank, ties.method = "min"), by = KEY]
setorderv(cb, c(KEY, "caller_rank"))
fwrite(cb, file.path(SC, "caller_rank.tsv"), sep = "\t")
print(dcast(cb, caller ~ gate + t, value.var = "caller_rank"))
print(cb[gate == "expected" & t == "t50", .(caller_rank, config, rank_sum, f1_arm, f1_sample, f1_celltype)])
