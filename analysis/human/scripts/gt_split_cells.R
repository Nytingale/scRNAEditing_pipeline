#!/usr/bin/env Rscript
# =============================================================================
# gt_split_cells.R   (2026-10-13)
# Split-by-cells circularity test (user 2026-10-13): NTC cells split at random into halves A and B within
# every teratoma x cell type (seed 20261013). Truth from half A (+ the KO arms, full), calls from half B,
# so truth-defining and called reads are disjoint.
# Input : _celltype_frags_umi_v6/_cells_NTC_ctrl.tsv (sample_id, cell_key, celltype; all 9,540 NTC cells),
#         _celltype_frags_umi_v6/_barcode_prefix_map.tsv (cell-type BAM CB prefix per teratoma)
# Output: gt_eval_v8/split_cells/barcodes/
#   cells_split.tsv                     sample_id, cell_key, celltype, half
#   T<n>_NTC_ctrl_<H>.tsv               barcode, celltype (gt_barcodes format; also the CB list of the per-teratoma BAMs)
#   merged_NTC_ctrl_<H>.tsv             the four teratomas pooled (raw barcodes)
#   NTC_ctrl__<celltype>_<H>.tsv        prefixed barcode, celltype (gt_barcodes_celltype format = CB of the cell-type BAMs)
#   split_counts.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
F <- file.path(PM, "_celltype_frags_umi_v6")
O <- file.path(PM, "gt_eval_v8", "split_cells", "barcodes")
dir.create(O, recursive = TRUE, showWarnings = FALSE)
SEED <- 20261013L

cells <- fread(file.path(F, "_cells_NTC_ctrl.tsv"))
pref <- fread(file.path(F, "_barcode_prefix_map.tsv"))
stopifnot(!anyDuplicated(cells[, .(sample_id, cell_key)]), all(cells$sample_id %in% pref$sample_id))
set.seed(SEED)
setorder(cells, sample_id, celltype, cell_key)
# random half within each teratoma x cell type; odd groups give the extra cell to A or B at random
cells[, half := { h <- rep(c("A", "B"), length.out = .N); if (.N %% 2 == 1 && runif(1) < 0.5) h <- chartr("AB", "BA", h); sample(h) },
      by = .(sample_id, celltype)]
cells[, prefixed := paste0(pref$code[match(sample_id, pref$sample_id)], "-", cell_key)]
fwrite(cells, file.path(O, "cells_split.tsv"), sep = "\t")

wr <- function(x, f) fwrite(x, file.path(O, f), sep = "\t", col.names = FALSE)
for (h in c("A", "B")) {
  for (s in unique(cells$sample_id)) wr(cells[sample_id == s & half == h, .(cell_key, celltype)], paste0(s, "_", h, ".tsv"))
  wr(cells[half == h, .(cell_key, celltype)], paste0("merged_NTC_ctrl_", h, ".tsv"))
  for (ct in unique(cells$celltype)) wr(cells[celltype == ct & half == h, .(prefixed, celltype)], sprintf("NTC_ctrl__%s_%s.tsv", ct, h))
}
cnt <- dcast(cells[, .N, by = .(sample_id, half)], sample_id ~ half, value.var = "N")
cnt <- rbind(cnt, dcast(cells[, .N, by = .(sample_id = paste0("NTC_ctrl__", celltype), half)], sample_id ~ half, value.var = "N"))
fwrite(cnt, file.path(O, "split_counts.tsv"), sep = "\t")
print(cnt); cat("cells A/B:", cells[, sum(half == "A")], "/", cells[, sum(half == "B")], "\n")
