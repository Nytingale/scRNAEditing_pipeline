#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_gtdefinition_source_data.R   (split from GTDEF_gtdefinition_plot.R 2026-10-25)
#
# Source data of the ground-truth threshold figures (GTDEF_step1_quality, GTDEF_step2_universe, GTDEF_robustness,
# GTDEF_robust_r1_centre_neighbours, GTDEF_robust_r5_informative_ko, GTDEF_panel_sweeps,
# GTDEF_supp_gt_thresholds_v8, GTDEF_v8_quality_universe). Each part is one report table, named without the
# GTDEF_ prefix:
#   step1_spectrum, step1_ntc_vs_ko, step1_summary      BQ/MQ 30 vs 20 substitution spectrum (every aligned base)
#   step2_summary, step2_af_hist, step2_ko_response     dbSNP common vs other sites (every candidate)
#   robust_centre_neighbours, robust_m_theory, robust_ko_cov_needed, robust_m_curve   one threshold at a time (v6)
#   robust_centre_neighbours_v8, robust_cleaning_v8, robust_m_curve_v8, robust_ko_cov_curve_v8   the same, v8 sites
#   sweep_caller_best{,_v8}, sweep_rank_stability{,_v8}  best configuration per caller vs depletion threshold
#
# INPUT : report/pm_gt/GTDEF_{step1,step2,robust}_*.tsv, GTDEF_sweep_*_v8.tsv   <- GTDEF_prep_source_data.R
#         human/pm/gt_eval_v6/scoring/sweep_infneg/{caller_best,sweep_rank_stability}.tsv  <- gt_sweep_infneg_v6.R
#           (copied to report/pm_gt/GTDEF_sweep_{caller_best,rank_stability}.tsv)
# OUTPUT: source_data/pm_gt/GTDEF_gtdefinition.<part>.tsv (drawn by GTDEF_gtdefinition_plot.R)
# ============================================================================
suppressPackageStartupMessages(library(data.table))
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_gtdefinition"
rp <- function(stem) file.path(REPORT_DIR, paste0(stem, ".tsv"))

SW <- file.path(GT_V6_DIR, "scoring", "sweep_infneg")
file.copy(file.path(SW, "caller_best.tsv"), rp("GTDEF_sweep_caller_best"), overwrite = TRUE)
file.copy(file.path(SW, "sweep_rank_stability.tsv"), rp("GTDEF_sweep_rank_stability"), overwrite = TRUE)

PARTS <- c("step1_spectrum", "step1_ntc_vs_ko", "step1_summary",
           "step2_summary", "step2_af_hist", "step2_ko_response",
           "robust_centre_neighbours", "robust_m_theory", "robust_ko_cov_needed", "robust_m_curve",
           "robust_centre_neighbours_v8", "robust_cleaning_v8", "robust_m_curve_v8", "robust_ko_cov_curve_v8",
           "sweep_caller_best", "sweep_rank_stability", "sweep_caller_best_v8", "sweep_rank_stability_v8")
for (part in PARTS) {
  f <- rp(paste0("GTDEF_", part))
  if (!file.exists(f)) { message("missing ", f); next }
  write_sd(fread(f, sep = "\t", showProgress = FALSE, integer64 = "double"), SRC_DIR, STEM, part)
}
message("source data -> ", SRC_DIR)
