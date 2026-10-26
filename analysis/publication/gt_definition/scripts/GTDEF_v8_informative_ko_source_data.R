#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_informative_ko_source_data.R   (split from GTDEF_gtdefinition_plot.R 2026-10-25)
# Source data for response figure 1 (the informative KO arm: m = 3 and the coverage gates). Panels a (e^-m) and
# c (KO coverage each gate requires) are theory, computed in the plot script.
#   b_null_fit    sites with 0 variant reads in the KO vs m: non-AG/TC and A-to-G candidate sites, fit and 95% band
#   d_m_curve     power-gated: sites passing and ground-truth sites at t = 0/25/50/75% vs m
#   e_cov_curve   depth-gated: the same vs KO coverage
# INPUT : report/pm_gt/GTDEF_v8_m_null_fit.tsv (GTDEF_v8_m_null_check_plot.R),
#         report/pm_gt/GTDEF_robust_{m_curve,ko_cov_curve}_v8.tsv (GTDEF_prep_source_data.R, from
#         human/pm/gt_eval_v8/robustness/)
# OUTPUT: source_data/pm_gt/GTDEF_v8_informative_ko.{b_null_fit,d_m_curve,e_cov_curve}.tsv
#         (drawn by GTDEF_v8_informative_ko_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_v8_informative_ko"
rd <- function(stem) {
  f <- file.path(REPORT_DIR, paste0(stem, ".tsv"))
  stopifnot(file.exists(f))
  fread(f, sep = "\t", showProgress = FALSE)
}

write_sd(rd("GTDEF_v8_m_null_fit"), SRC_DIR, STEM, "b_null_fit")
write_sd(rd("GTDEF_robust_m_curve_v8"), SRC_DIR, STEM, "d_m_curve")
write_sd(rd("GTDEF_robust_ko_cov_curve_v8"), SRC_DIR, STEM, "e_cov_curve")
