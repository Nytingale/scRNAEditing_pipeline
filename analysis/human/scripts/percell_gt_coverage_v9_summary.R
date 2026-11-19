#!/usr/bin/env Rscript
# percell_gt_coverage_v9_summary.R   (v9 copy 2026-11-17 of the v8 summary, 2026-10-15; review R1.2 / R2.2 per-cell sparsity, v9 truth)
# Per-cell depth at the v9 KO ground-truth positives, same quantities as PERCELL_gt_site_coverage_summary.tsv (v5):
# cell x site pairs with >= 1 read and >= MIN_COV reads; per site x sample, the best-covered single cell and the
# sample pseudobulk (all its cells) reaching MIN_COV. One row per truth (gate x t) and for the union of positives.
# Input : human/pm/percell_gt_coverage_v9/{<S>_gt_site_counts.tsv.gz, gt_v9_site_classes.tsv}
#         (human/pm/scripts/percell_gt_site_counts_v9.sh; unprocessed reads, NH == 1, BQ / MQ >= 30, no dedup)
# Output: human/pm/percell_gt_coverage_v9/percell_gt_coverage_v9_summary.tsv
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
D <- file.path(PM, "percell_gt_coverage_v9")
MIN_COV <- 10L

files <- list.files(D, pattern = "_gt_site_counts\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(files) == 20)
cnt <- rbindlist(lapply(files, function(f) fread(cmd = paste("zcat", shQuote(f)))))
cnt[, depth := ref_reads + alt_reads]
cells <- fread(file.path(PM, "seurat_analysis", "ground_truth_cells_paper.tsv"))[, .(sample = paste0(sample, "_", arm), barcode)]
cnt <- cnt[cells, on = .(sample, barcode), nomatch = 0L]
N_CELLS <- nrow(cells); SAMPLES <- unique(cells$sample)

cls <- melt(fread(file.path(D, "gt_v9_site_classes.tsv")), id.vars = c("gate", "site"), variable.name = "t", value.name = "class")
cls <- cls[class == "positive"][, t := sub("class_", "", t)]
sets <- c(split(cls$site, paste(cls$gate, cls$t)), list(`union` = unique(cls$site)))

summ <- rbindlist(lapply(names(sets), function(nm) {
  s <- sets[[nm]]; x <- cnt[site %in% s]; n_pairs <- as.numeric(N_CELLS) * length(s)
  ps <- x[, .(best = max(depth), pb = sum(depth)), by = .(sample, site)]
  n_ss <- length(SAMPLES) * length(s)
  data.table(truth = nm, cells = N_CELLS, sites = length(s), cell_site_pairs = n_pairs,
             pct_pairs_ge1_read = 100 * nrow(x) / n_pairs,
             pct_pairs_ge_min_cov = 100 * x[depth >= MIN_COV, .N] / n_pairs,
             pct_site_sample_best_cell_ge_min_cov = 100 * ps[best >= MIN_COV, .N] / n_ss,
             pct_site_sample_pseudobulk_ge_min_cov = 100 * ps[pb >= MIN_COV, .N] / n_ss)
}))
print(summ, digits = 4)
fwrite(summ, file.path(D, "percell_gt_coverage_v9_summary.tsv"), sep = "\t")
