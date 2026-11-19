#!/usr/bin/env Rscript
# gt_v8_celltype_covered_positives.R   (2026-10-30)
# Covered positives (n_pos_unit: >= 10 stranded reads in the unit) per cell type under the nine v8 truths, from a
# mismatch-state scoring of all 23 NTC cell types (run_gt_v8_rescore.sh task 6 and 14 with GT_UNIT_KEEP= empty).
# The input of gt_v8_celltype_keep_pw.R (written inline on 2026-10-13; same code).
# Usage : Rscript gt_v8_celltype_covered_positives.R <KO scoring dir> <REDIportal scoring dir>
#         e.g. gt_eval_v8/_ct_all23/scoring/mismatch gt_eval_v8/_ct_all23/rediportal/scoring_mismatch
# Output: gt_eval_v8/unit_cov/celltype_covered_positives.tsv
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 2)
V8 <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
k <- unique(fread(file.path(a[1], "celltype", "gt_v7_unit_scores.tsv"))[, .(truth = paste(gate, t, sep = "_"), unit, n = n_pos_unit)])
r <- unique(fread(file.path(a[2], "unit_scores.tsv"))[tier == "celltype", .(truth = paste0("redi_", variant), unit, n = n_pos_unit)])
x <- rbind(k, r); x[, unit := sub("NTC_ctrl__", "", unit)]
w <- dcast(x, unit ~ truth, value.var = "n")
setcolorder(w, c("unit", "expected_t25", "expected_t50", "expected_t75", "depth_t25", "depth_t50", "depth_t75",
                 "redi_cov10_alt1", "redi_alt3", "redi_alt3_e05"))
stopifnot(nrow(w) == 23, !anyNA(w))
w[, min_all := do.call(pmin, .SD), .SDcols = -1]; setorder(w, min_all)
print(w)
fwrite(w, file.path(V8, "unit_cov", "celltype_covered_positives.tsv"), sep = "\t")
