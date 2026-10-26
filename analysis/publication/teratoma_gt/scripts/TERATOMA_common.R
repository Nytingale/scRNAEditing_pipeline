#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_common.R
# Shared paths + figure-saving wrapper for the teratoma ground-truth RESULTS
# deliverable. Sourced by TERATOMA_prep_source_data.R and TERATOMA_results_plot.R.
# Do not put figure code here.
# ============================================================================
# SCOPE. This deliverable holds the BENCHMARK RESULT: how each method
# combination scores against the ground truth. The panels that justify the
# ground-truth DEFINITION live separately in publication/gt_definition/ — a
# definition and a result are different claims and should not share a folder.
#
# The five `FigS_teratoma_*` panels here are marked `naming = "legacy"` in
# FIGURE_INDEX.tsv and are superseded by the `TERATOMA_*` panels; they are left on
# disk rather than deleted, as nothing in this project is removed silently.
#
# Naming convention (see ../../README.md):
#   TERATOMA_<plottype>_<description>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

teratoma_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) return(dirname(normalizePath(f[1], mustWork = FALSE)))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(o)) return(dirname(normalizePath(o, mustWork = FALSE)))
  getwd()
}

SCRIPT_DIR <- teratoma_script_dir()
DELIV_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_DIR   <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

source(file.path(REPO_DIR, "plotting", "plotting_functions.R"))
# v8 paper truth (2026-10-13): V8_* constants, human/pm/gt_eval_v8/README.md
source(file.path(REPO_DIR, "plotting", "gt_v8_paths.R"))

# Per-dataset output tree, so a second dataset (Replogle) cannot overwrite the
# first's tables and figures. See publication/README.md.
teratoma_set_dataset <- function(ds = "pm_gt") {
  DATASET    <<- ds
  FIG_DIR    <<- file.path(DELIV_DIR, "figures",     ds)
  SRC_DIR    <<- file.path(DELIV_DIR, "source_data", ds)
  REPORT_DIR <<- file.path(DELIV_DIR, "report",      ds)
  for (d in c(FIG_DIR, SRC_DIR, REPORT_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(ds)
}

# ── Inputs ─────────────────────────────────────────────────────────────────
# Current truth (2026-09-29): the v6 sweep scoring, positives at t = 25/50/75%
# and fixed informative negatives (human/pm/gt_eval_v6/GROUND_TRUTH_DEFINITION.md),
# written in the v5 table layout by human/scripts/gt_teratoma_plot_inputs_v6.R.
# ag_percent is already in the scoring table (from the scored call files).
# The v5 inputs are the fallback only; pre-v6 copy of this file:
# ../_superseded/v5_results_2026-09-29/scripts/TERATOMA_common.R
GT_DIR      <- file.path(REPO_DIR, "human", "pm", "gt_eval")
SCORING_DIR <- file.path(REPO_DIR, "human", "pm", "gt_eval_v6", "scoring", "sweep_infneg", "teratoma_inputs")
if (!file.exists(file.path(SCORING_DIR, "gt_tier_scoring.tsv"))) {
  SCORING_DIR <- file.path(GT_DIR, "true_def_v5_full")
  message("WARNING: v6 sweep scoring not found; using v5 (true_def_v5_full)")
}
if (!file.exists(file.path(SCORING_DIR, "gt_tier_scoring.tsv")))
  stop("missing ", file.path(SCORING_DIR, "gt_tier_scoring.tsv"))
message("scoring dir: ", SCORING_DIR)

# The same tables for the three mismatch-window-filter states (HANDOFF 2026-09-29 UPDATE-22):
# none = step 4 (no window filter), mismatch = step 5 (nearby-mismatch filter only),
# both = step 6 (the benchmark call set, SCORING_DIR).
WINDOW_SCORING_DIRS <- c(
  none     = file.path(REPO_DIR, "human", "pm", "gt_eval_v6", "prewindow", "scoring", "sweep_infneg", "teratoma_inputs"),
  mismatch = file.path(REPO_DIR, "human", "pm", "gt_eval_v6", "prewindow", "scoring_step5", "sweep_infneg", "teratoma_inputs"),
  both     = SCORING_DIR)
# States shown in the cross-tier grids: the full benchmark call set only (HANDOFF UPDATE-25).
# Without window filters every call set is >= 31% non-A-to-G; with step 5 only the median
# is 58-70% A-to-G (TERATOMA_mismatch_composition_plot.R, HANDOFF UPDATE-24). The Local mismatch
# context column is drawn only when more than one state is shown.
GRID_WINDOW_STATES <- "both"
GRID_SHOW_WINDOW <- length(GRID_WINDOW_STATES) > 1
read_window_scoring <- function(rel, states = names(WINDOW_SCORING_DIRS)) {
  rbindlist(lapply(states, function(w) {
    f <- file.path(WINDOW_SCORING_DIRS[[w]], rel)
    if (!file.exists(f)) stop("missing ", f)
    fread(f)[, window_filter := w]
  }))
}

# Nine global truth sets used for ranking robustness (2026-10-01). Ground truth: NTC cov >= 10,
# alt >= 3, edit >= 0.05, not dbSNP common, reproducible across teratomas, and a covered KO arm:
#   expected  KO cov x NTC edit >= 3 (frozen v6; sweep_infneg)
#   depth     NTC and KO cov >= 10 (sweep_infneg_depth)
# REDIportal (not dbSNP common, merged NTC): cov >= 10 with alt >= 1 / alt >= 3 / alt >= 3 and edit >= 5%.
TRUTH_DEPTH_DIR <- file.path(REPO_DIR, "human", "pm", "gt_eval_v6", "scoring", "sweep_infneg_depth", "teratoma_inputs")
TRUTH_REDI_DIR  <- file.path(REPO_DIR, "human", "pm", "gt_eval_v6", "rediportal", "scoring")
NINE_TRUTHS <- data.table(
  truth = c(paste0("gt_expected_", c("d0.25", "d0.50", "d0.75")), paste0("gt_depth_", c("d0.25", "d0.50", "d0.75")),
            paste0("redi_", c("cov10_alt1", "alt3", "alt3_e05"))),
  label = c(sprintf("Ground truth, \u2265 %d%% depletion, expected KO alt \u2265 3", c(25, 50, 75)),
            sprintf("Ground truth, \u2265 %d%% depletion, NTC and KO coverage \u2265 10", c(25, 50, 75)),
            "REDIportal, \u2265 10 reads, \u2265 1 G", "REDIportal, \u2265 10 reads, \u2265 3 G",
            "REDIportal, \u2265 10 reads, \u2265 3 G, \u2265 5% editing level"),
  colour = unname(V8_TRUTH_COLOURS),
  file = c(file.path(SCORING_DIR, "invariance", c("d0.25", "d0.50", "d0.75"), "gt_tier_scoring.tsv"),
           file.path(TRUTH_DEPTH_DIR, "invariance", c("d0.25", "d0.50", "d0.75"), "gt_tier_scoring.tsv"),
           file.path(TRUTH_REDI_DIR, paste0("gt_tier_scoring_", c("cov10_alt1", "alt3", "alt3_e05"), ".tsv"))))
read_truth_scoring <- function(id) {
  f <- NINE_TRUTHS[truth == id, file]
  if (!length(f) || !file.exists(f)) stop("missing truth scoring for ", id, ": ", f)
  fread(f)
}

SCORING_TSV    <- file.path(SCORING_DIR, "gt_tier_scoring.tsv")
SCORING_SE_TSV <- file.path(SCORING_DIR, "gt_tier_scoring_with_se.tsv")
DISC_SAMPLE_TSV <- file.path(GT_DIR, "persample_union_v4", "gt_eval_discovery.tsv.gz")
COLLATED       <- file.path(REPO_DIR, "human", "pm", "results_gt_persample",
                            "summary", "collated", "all_results_summary.tsv")

# The grid dimensions that define a unique method combination. Reporting by
# caller alone is not answerable: a reviewer immediately asks which dedup,
# which stranding, which SNP treatment, and each gives different numbers.
COMBO_DIMS <- c("caller", "bam_type", "strand_mode", "snp_subtract")

save_figure <- function(plot, stem, width, height, sources, data = NULL,
                        grid = NA_character_, rasterized = FALSE,
                        script = NA_character_) {
  svg <- save_panel(plot, file.path(FIG_DIR, stem),
                    width = width, height = height,
                    sources = sources, data = data, grid = grid,
                    rasterized = rasterized, script = script,
                    sidecar_dir = SRC_DIR)
  message("  figure      : ", file.path("figures", DATASET, basename(svg)))
  invisible(svg)
}
