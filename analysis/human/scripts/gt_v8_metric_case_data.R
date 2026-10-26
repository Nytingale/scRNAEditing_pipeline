#!/usr/bin/env Rscript
# ============================================================================
# gt_v8_metric_case_data.R   (2026-10-22)
# Data of Response Fig 8 (TERATOMA_v8_metric_case_plot.R), replacing gt_v8_sites_inflation.R (old filter sweep and
# old scoring, REDItools3 without its 10-read floor). Per configuration (58) and setting:
#   centre  published step-5 call sets (results_gt_v8, REDItools3 floor applied) + headline scoring
#           (gt_eval_v8/scoring_pwcov/mismatch/config_rank.tsv)
#   alt1    filter sweep 2, >= 3 variant-read minimum removed (filter_sweep2/results/alt1) + the same scoring
#           (filter_sweep2/scoring_pwcov_extra/alt1/mismatch/config_rank.tsv)
# A-to-G sites / percentage = merged NTC arm (NTC_ctrl) calls; mean F1 / precision / recall = mean over the 6 KO-based
# truths of the cross-tier means.
# Output: human/pm/gt_eval_v8/metric_case/metric_case.tsv
# ============================================================================
suppressPackageStartupMessages(library(data.table))
V8 <- "human/pm/gt_eval_v8"
OUT <- file.path(V8, "metric_case"); dir.create(OUT, showWarnings = FALSE)
SRC <- list(
  centre = list(calls = "human/pm/results_gt_v8/NTC_ctrl", rank = file.path(V8, "scoring_pwcov/mismatch/config_rank.tsv")),
  alt1 = list(calls = file.path(V8, "filter_sweep2/results/alt1/results_gt_v8/NTC_ctrl"),
              rank = file.path(V8, "filter_sweep2/scoring_pwcov_extra/alt1/mismatch/config_rank.tsv")))
cfg_of <- function(f) {
  p <- strsplit(f, "/")[[1]]; i <- which(p == "04_filtered")
  paste(p[i + 1], p[i + 2], p[i + 3], sub("snp_", "", p[i + 4]), sep = "__")
}
d <- rbindlist(lapply(names(SRC), function(s) {
  r <- fread(SRC[[s]]$rank)
  stopifnot(r[, uniqueN(paste(gate, t))] == 6L)
  sc <- r[, .(mean_f1 = mean(mean_f1),
              precision = mean((precision_arm + precision_sample + precision_celltype) / 3),
              recall = mean((recall_arm + recall_sample + recall_celltype) / 3)), by = .(config, caller)]
  fs <- list.files(SRC[[s]]$calls, "_step5_nomm\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  cnt <- rbindlist(lapply(fs, function(f) {
    v <- fread(cmd = paste("zcat", shQuote(f)), select = "varType")$varType
    data.table(config = cfg_of(f), n = length(v), ag = sum(v == "AG"))
  }))
  x <- merge(sc, cnt, by = "config")
  stopifnot(nrow(x) == 58L)
  x[, `:=`(setting = s, pct = 100 * ag / n)]
}))
setcolorder(d, c("setting", "config", "caller", "n", "ag", "pct"))
fwrite(d, file.path(OUT, "metric_case.tsv"), sep = "\t")
print(d[, .(configs = .N, median_ag = median(ag), median_pct = median(pct), median_f1 = median(mean_f1),
            rho_sites = cor(ag, mean_f1, method = "spearman"), rho_pct = cor(pct, mean_f1, method = "spearman")),
        by = setting])
