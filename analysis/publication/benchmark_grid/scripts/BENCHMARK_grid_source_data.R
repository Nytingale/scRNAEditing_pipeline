#!/usr/bin/env Rscript
# ============================================================================
# BENCHMARK_grid_source_data.R   (split from BENCHMARK_grid_plot.R 2026-10-25)
#
# Source data of the benchmark-grid panels (caller x bam_type x strand_mode x snp_subtract, per sample):
#   a_final          final A-to-G percentage per run                       (BENCHMARK_box_ag_pct_final)
#   b_attrition      surviving sites per filtering step, median and IQR    (BENCHMARK_line_site_attrition)
#   c_ag_pct         A-to-G percentage per filtering step, median and IQR  (BENCHMARK_line_ag_pct_progression)
#   d_strand_shift   hybrid minus annotation A-to-G percentage, paired within run  (BENCHMARK_point_strand_mode_shift)
#   e_windows        adaptive filter windows and step-6 A-to-G percentage  (BENCHMARK_heatmap_window)
#   f_strand_source  final sites by strand source in hybrid mode           (BENCHMARK_stackedbar_strand_source)
# d-f are absent for datasets run with one strand mode. BENCHMARK_panel_grid is drawn from a_final and c_ag_pct.
# Every part keeps each grid dimension; assert_no_dropped_combos() checks that only `sample` is aggregated.
#
# INPUT : data/<DATASET>/BENCHMARK_{final_by_run,attrition_by_step}.tsv.gz,
#         BENCHMARK_{strand_mode_cmp,windows,strand_source}.tsv   <- BENCHMARK_prep_source_data.R
# OUTPUT: source_data/<DATASET>/BENCHMARK_grid.<part>.tsv (drawn by BENCHMARK_grid_plot.R)
# Usage : Rscript BENCHMARK_grid_source_data.R [DATASET]   (default pm_gt_persample)
# ============================================================================
suppressPackageStartupMessages(library(data.table))
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(SCRIPT_DIR, "BENCHMARK_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
args <- commandArgs(trailingOnly = TRUE)
bench_set_dataset(if (length(args) >= 1) args[1] else "pm_gt_persample")
STEM <- "BENCHMARK_grid"

need <- function(p) {
  if (!file.exists(p)) stop("missing prepared input: ", p, "\nRun: Rscript scripts/BENCHMARK_prep_source_data.R")
  p
}
for (part in c("d_strand_shift", "e_windows", "f_strand_source"))
  unlink(c(sd_file(SRC_DIR, STEM, part), sd_schema_file(SRC_DIR, STEM, part)))

dat.a <- fread(need(file.path(DATA_DIR, "BENCHMARK_final_by_run.tsv.gz")))
dat.a[, bam_type := factor_keep_all(bam_type, bam_type_levels, "bam_type")]
assert_no_dropped_combos(dat.a, shown = GRID_DIMS, panel = "BENCHMARK_box_ag_pct_final", aggregated = AGG_DIMS)
write_sd(dat.a, SRC_DIR, STEM, "a_final")

attr_all <- fread(need(file.path(DATA_DIR, "BENCHMARK_attrition_by_step.tsv.gz")))
attr_all[, step_label := factor(step_label, levels = STEP_LABELS)]
attr_all[, bam_type := factor_keep_all(bam_type, bam_type_levels, "bam_type")]
for (x in list(c("b_attrition", "n_sites", "BENCHMARK_line_site_attrition"),
               c("c_ag_pct", "ag_pct", "BENCHMARK_line_ag_pct_progression"))) {
  d <- attr_all[metric == x[2]]
  assert_no_dropped_combos(d, shown = c(GRID_DIMS, "step_num"), panel = x[3], aggregated = AGG_DIMS,
                           position = c(GRID_DIMS, "step_num"))
  write_sd(d, SRC_DIR, STEM, x[1])
}

f <- file.path(DATA_DIR, "BENCHMARK_strand_mode_cmp.tsv")
if (file.exists(f)) {
  dat.d <- fread(f)
  dat.d[, bam_type := factor_keep_all(bam_type, bam_type_levels, "bam_type")]
  assert_no_dropped_combos(dat.d, shown = setdiff(GRID_DIMS, "strand_mode"), panel = "BENCHMARK_point_strand_mode_shift",
                           aggregated = c(AGG_DIMS, "strand_mode"))
  write_sd(dat.d, SRC_DIR, STEM, "d_strand_shift")
}

f <- file.path(DATA_DIR, "BENCHMARK_windows.tsv")
if (file.exists(f)) {
  dat.e <- fread(f)
  dat.e[, bam_type := factor_keep_all(bam_type, bam_type_levels, "bam_type")]
  assert_no_dropped_combos(dat.e, shown = GRID_DIMS, panel = "BENCHMARK_heatmap_window", aggregated = AGG_DIMS,
                           position = GRID_DIMS)
  write_sd(dat.e, SRC_DIR, STEM, "e_windows")
}

f <- file.path(DATA_DIR, "BENCHMARK_strand_source.tsv")
if (file.exists(f)) {
  dat.f <- fread(f)
  dat.f[, bam_type := factor_keep_all(bam_type, bam_type_levels, "bam_type")]
  dat.f[, source := factor(source, levels = c("read", "annot"))]
  assert_no_dropped_combos(dat.f, shown = c(setdiff(GRID_DIMS, "strand_mode"), "source"),
                           panel = "BENCHMARK_stackedbar_strand_source", aggregated = c(AGG_DIMS, "strand_mode"),
                           position = c(setdiff(GRID_DIMS, "strand_mode"), "source"))
  write_sd(dat.f, SRC_DIR, STEM, "f_strand_source")
}
message("source data -> ", SRC_DIR)
