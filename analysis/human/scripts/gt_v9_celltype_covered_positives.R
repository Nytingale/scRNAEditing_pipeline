#!/usr/bin/env Rscript
# gt_v9_celltype_covered_positives.R   (2026-11-16)
# v9 counterpart of gt_v8_celltype_covered_positives.R + gt_v8_celltype_keep_pw.R (gate_cov): covered positives
# (n_pos_unit, the per-unit recall denominator of the headline scoring) of each of the 23 NTC cell types under the
# nine v9 truths, from scorings of all 23 cell types (gt_v9_celltype_keep.sh STAGE=ko / redi, GT_UNIT_KEEP empty),
# and the cell types with >= MIN under every truth.
# Usage : Rscript gt_v9_celltype_covered_positives.R [MIN (30)]
# Input : gt_eval_v9/celltype_all23/ko/celltype/gt_v7_unit_scores.tsv, gt_eval_v9/celltype_all23/redi/unit_scores.tsv
# Output: gt_eval_v9/unit_cov/celltype_covered_positives_pwcov.tsv, celltype_keep_pwcov.txt;
#         compare with gt_eval_v8/unit_cov/celltype_keep_pwcov.txt (the list used by the v9 scorings) printed
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
MIN <- if (length(a)) as.integer(a[1]) else 30L
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
D <- file.path(PM, "gt_eval_v9", "celltype_all23")
k <- unique(fread(file.path(D, "ko", "celltype", "gt_v7_unit_scores.tsv"))[, .(truth = paste(gate, t, sep = "_"), unit, n = n_pos_unit)])
r <- unique(fread(file.path(D, "redi", "unit_scores.tsv"))[tier == "celltype", .(truth = paste0("redi_", variant), unit, n = n_pos_unit)])
x <- rbind(k, r); x[, unit := sub("NTC_ctrl__", "", unit)]
stopifnot(!anyDuplicated(x[, .(truth, unit)]))
w <- dcast(x, unit ~ truth, value.var = "n")
setcolorder(w, c("unit", "expected_t25", "expected_t50", "expected_t75", "depth_t25", "depth_t50", "depth_t75",
                 "redi_cov10_alt1", "redi_alt3", "redi_alt3_e05"))
stopifnot(nrow(w) == 23, !anyNA(w))
w[, min_all := do.call(pmin, .SD), .SDcols = -1]; setorder(w, -min_all, unit)
print(w)
keep <- paste0("NTC_ctrl__", w[min_all >= MIN, unit])
old <- readLines(file.path(PM, "gt_eval_v8", "unit_cov", "celltype_keep_pwcov.txt"))
cat(sprintf("v9: %d cell types with >= %d covered positives under every truth; v8 list %d; same set: %s\n",
            length(keep), MIN, length(old), setequal(keep, old)))
if (!setequal(keep, old)) cat("only v9:", setdiff(keep, old), "\nonly v8:", setdiff(old, keep), "\n")
dir.create(file.path(PM, "gt_eval_v9", "unit_cov"), showWarnings = FALSE)
fwrite(w, file.path(PM, "gt_eval_v9", "unit_cov", "celltype_covered_positives_pwcov.tsv"), sep = "\t")
writeLines(keep, file.path(PM, "gt_eval_v9", "unit_cov", "celltype_keep_pwcov.txt"))
