#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_window_state_rank_scatter_source_data.R   (2026-10-18; split from the _plot.R script 2026-10-25)
# Source data for response figure 7: does the removal of sites near other mismatch types change the
# configuration ranking? Overall cross-tier rank (gt_rank_v7.R) of the 58 ranked configurations without the
# local context filter (state none = step 4) against the main analysis (state mismatch = step 5), t = 50% KO truths.
# Input : human/pm/gt_eval_v8/state_summary/config_ranks.tsv (human/scripts/gt_v8_state_caller_summary.R)
# Output: source_data/pm_gt/TERATOMA_v8_window_state_rank_scatter.{ranks,top_k,rho}.tsv
#         (drawn by TERATOMA_v8_window_state_rank_scatter_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_window_state_rank_scatter"

SRC <- c(rank = file.path(V8_DIR, paste0("state_summary", V8_SFX), "config_ranks.tsv"))
r <- fread(SRC[["rank"]])[state %in% c("none", "mismatch")]
stopifnot(setequal(r$truth, V8_TRUTH_IDS), r[, .N, by = .(state, truth)][, all(N == 58L)])
w <- dcast(r, config + caller + bam_type + truth ~ state, value.var = "overall_rank")
print(w[, .(rho = round(cor(none, mismatch, method = "spearman"), 3)), by = truth])
PLOT_TRUTHS <- c("expected_t50", "depth_t50")
w <- w[truth %in% PLOT_TRUTHS]
TRUTH_STRIP <- setNames(gsub(", ", ",\n", V8_TRUTH_LABELS[PLOT_TRUTHS], fixed = TRUE), PLOT_TRUTHS)
w[, `:=`(Truth = factor(TRUTH_STRIP[truth], TRUTH_STRIP), method = caller_factor(caller),
         proc = relabel(bam_type, label_bamType))]
PROC_LEVELS <- intersect(label_bamType[c("sampleRmDup", "sampleUmiDup", "recal", "recalq20", "recaldefault")], w$proc)
w[, proc := factor(proc, PROC_LEVELS)]
st <- w[, .(rho = cor(none, mismatch, method = "spearman")), by = Truth]
print(st[, .(Truth = gsub("\n", " ", Truth), rho = round(rho, 3))])
# top box: smallest k (>= 3) at which the top-k configurations are the same set with and without the filter
top_k <- w[, .(k = Filter(function(k) setequal(config[none <= k], config[mismatch <= k]), 3:58)[1]), by = Truth]
print(top_k)
top_k[, lab := sprintf("Same top %d\nconfigurations", k)]

write_sd(w, SRC_DIR, STEM, "ranks")
write_sd(top_k, SRC_DIR, STEM, "top_k")
write_sd(st, SRC_DIR, STEM, "rho")
