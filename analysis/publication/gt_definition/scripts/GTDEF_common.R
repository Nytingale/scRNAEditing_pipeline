#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_common.R
# Shared paths + figure-saving wrapper for the ground-truth-definition
# deliverable. Sourced by GTDEF_prep_source_data.R and the GTDEF_*_source_data.R scripts.
# Do not put figure code here.
# ============================================================================
# WHAT THIS DELIVERABLE IS FOR
# Every threshold in the v3 ground-truth definition
# (human/pm/gt_eval/true_def_v3/GROUND_TRUTH_DEFINITION.md) needs a panel that
# shows the measurement behind it, so a reviewer can see the choice was made
# from the data rather than asserted. It also carries the two negative results
# that changed the definition:
#   * the >=3-of-4 consistency criterion is DEPTH-CONFOUNDED, which is why
#     gt_eval_classify_v3.R defaults to --consistency_mode evaluable;
#   * --alt_rule binom is EXACTLY equivalent to count with min_alt = 2 at these
#     depths, which is the argument for min_alt = 2 over 3 and for keeping the
#     binomial machinery available but out of the benchmark grid.
#
# Naming convention (see ../../README.md):
#   GTDEF_<plottype>_<description>   one stem shared by the figure and its
#                                    source-data TSV.
# GTDEF_gtdefinition_source_data.R gathers the threshold tables and
# GTDEF_gtdefinition_plot.R draws every lettered panel from them.
#
# figures/ holds SVG only. save_figure() routes every TSV save_panel() emits
# into source_data/ instead.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

gtdef_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) return(dirname(normalizePath(f[1], mustWork = FALSE)))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(o)) return(dirname(normalizePath(o, mustWork = FALSE)))
  getwd()
}

SCRIPT_DIR <- gtdef_script_dir()
DELIV_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_DIR   <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

source(file.path(REPO_DIR, "plotting", "plotting_functions.R"))
# v8 paper truth (2026-10-13): V8_* constants, human/pm/gt_eval_v8/README.md
source(file.path(REPO_DIR, "plotting", "gt_v8_paths.R"))

# ── Per-dataset output tree ────────────────────────────────────────────────
# Only teratoma can support this analysis today (Replogle is K562: one cell
# line, three arms, no replicates and no cell types, so it has no tier ladder
# and no per-teratoma consistency to measure). The per-dataset layout is used
# anyway so a second dataset can never silently overwrite the first — the same
# latent bug benchmark_grid had. See publication/README.md.
gtdef_set_dataset <- function(ds = "pm_gt") {
  DATASET    <<- ds
  FIG_DIR    <<- file.path(DELIV_DIR, "figures",     ds)
  SRC_DIR    <<- file.path(DELIV_DIR, "source_data", ds)
  DATA_DIR   <<- file.path(DELIV_DIR, "data",        ds)
  REPORT_DIR <<- file.path(DELIV_DIR, "report",      ds)
  for (d in c(FIG_DIR, SRC_DIR, DATA_DIR, REPORT_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(ds)
}

# ── Inputs, resolved from the repo root ────────────────────────────────────
GT_DIR        <- file.path(REPO_DIR, "human", "pm", "gt_eval")
TRUTH_PILEUPS <- file.path(GT_DIR, "truth_combined")
UNION_DIR     <- file.path(GT_DIR, "persample_union")
TIER_DEPTH    <- file.path(GT_DIR, "tier_depth.tsv")
TRUE_DEF_V3   <- file.path(GT_DIR, "true_def_v3")
STRAND_TEST   <- file.path(REPO_DIR, "nextflow", "rendercheck", "caller_strand_test")

# Thresholds under examination. Held here so the prep and plot stages cannot
# disagree about which line is being drawn, and so changing a default is a
# one-line edit that propagates to every panel and caption.
GT_MIN_COV        <- 10
GT_MIN_EDIT_READS <- 3
GT_MIN_CTRL_EDIT  <- 0.05
GT_CONSISTENCY    <- 3
GT_TERATOMAS      <- c("T1", "T2", "T3", "T4")
GT_CTRL_ARM       <- "NTC_ctrl"
GT_TREAT_ARMS     <- c("ADAR1_KO", "ADARB1_KO")   # ADARB2/ADAR3 excluded, §2a.5
GT_DEPLETION      <- 0.50
# The four levels §2a.4 commits to reporting. Used where depletion is a
# CATEGORICAL aesthetic (the power panel colours), because nine hues would be
# unreadable.
GT_DEP_SERIES     <- c(0.25, 0.50, 0.75, 0.90)
# The finer grid used where depletion is a CONTINUOUS axis (the depletion
# series panel), so the curve is evenly sampled at 10% intervals.
GT_DEP_BREAKS     <- seq(0.1, 0.9, 0.1)
GT_SEQ_ERROR      <- 0.001
GT_ALPHA          <- 0.05

# ── GT v6 (genome-wide, caller-independent) ────────────────────────────────
# Used by GTDEF_prep_source_data.R / GTDEF_gtdefinition_source_data.R (v6 versions).
# The v3/v5 constants above are kept because GTDEF_aavs1_specificity_plot.R and
# GTDEF_bam_processing_plot.R still read v5 inputs through this file.
# Pre-edit copy of this file: ../_superseded/scripts/GTDEF_common.R
GT_V6_DIR     = file.path(REPO_DIR, "human", "pm", "gt_eval_v6")
GT_V6_SWEEPS  = file.path(GT_V6_DIR, "sweeps")
GT_V6_SCAN    = file.path(GT_V6_DIR, "scan")
GT_V6_PILEUPS = file.path(GT_V6_DIR, "pileups")
GT_V6_BAMDIR  = file.path(REPO_DIR, "human", "pm", "02_gt_bams_umi")
GT_V6_FROZEN  = file.path(GT_V6_DIR, "FROZEN_RULES.tsv")
GT_V6_STEPS   = c("step1_quality", "step2_universe", "step3_coverage",
                  "step4_control_editing", "step5_reproducibility", "step6_depletion")

# Current v6 thresholds: FROZEN_RULES.tsv where frozen, else the provisional
# values written by run_gt_truth_rebuild_v6.sh (sweeps/provisional.tsv).
# The single source of parameter values for the v6 GTDEF scripts.
gt_v6_params = function() {
  prov = data.table::fread(file.path(GT_V6_SWEEPS, "provisional.tsv"), colClasses = "character")
  value = setNames(prov$value, prov$parameter)
  frozen = setNames(rep(FALSE, length(value)), names(value))
  if (file.exists(GT_V6_FROZEN)) {
    fr = data.table::fread(GT_V6_FROZEN, colClasses = "character")
    value[fr$parameter] = fr$value
    frozen[fr$parameter] = TRUE
  }
  list(value = value, frozen = frozen)
}
# "0.1" == "0.10", "evaluable" == "evaluable"
gt_same_value = function(a, b) {
  na = suppressWarnings(as.numeric(a)); nb = suppressWarnings(as.numeric(b))
  ifelse(!is.na(na) & !is.na(nb), abs(na - nb) < 1e-9, as.character(a) == as.character(b))
}

# ── Figure writer ──────────────────────────────────────────────────────────
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
