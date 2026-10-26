#!/usr/bin/env Rscript
# WASP_common.R
# Shared paths + figure-saving wrapper for the WASP reference-bias supplement.
# Sourced by the WASP_*_source_data.R scripts. Do not put figure code here.
#
# Naming convention for this deliverable (see ../README.md):
#   WASP_<plottype>_<description>   one stem shared by the figure and its
#                                   source-data TSV.
# WASP_refbias_plot.R builds every lettered panel and the combined figure they
# assemble into, from the tables of WASP_refbias_source_data.R. The WASP_ prefix groups them together in
# a listing; figure numbers are assigned later from publication/FIGURE_INDEX.tsv.
#
# EVERY PATH IS PER DATASET. The diagnostic runs on more than one dataset
# (Reviewer 2, Major 4 asks for a second), so data/, report/, figures/ and
# source_data/ each carry a <dataset> subfolder. Figure STEMS stay identical
# across datasets — the folder says which dataset, so panels stay directly
# comparable and no stem has to be renamed when a dataset is added.
#
# figures/<dataset>/ holds SVG only. save_figure() routes every TSV
# save_panel() emits (<stem>.sources.tsv, <stem>.source_data.tsv,
# figure_source_map.tsv) into source_data/<dataset>/ instead.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ── Locate the repo so this deliverable can be run from anywhere ────────────
wasp_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) return(dirname(normalizePath(f[1], mustWork = FALSE)))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(o)) return(dirname(normalizePath(o, mustWork = FALSE)))
  getwd()
}

SCRIPT_DIR <- wasp_script_dir()
DELIV_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_DIR   <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

source(file.path(REPO_DIR, "plotting", "plotting_functions.R"))

# ── Dataset registry ────────────────────────────────────────────────────────
# WASP_datasets.tsv is the single source of truth, shared with the bash stages
# so a dataset is described once. Adding a dataset is a row, not a code change.
WASP_DATASETS <- file.path(SCRIPT_DIR, "WASP_datasets.tsv")

wasp_dataset <- function(ds) {          # not `dataset`: that is a column name
  reg <- fread(WASP_DATASETS)
  row <- reg[dataset == ds]
  if (nrow(row) != 1L) {
    stop("Unknown dataset '", ds, "'. Known: ",
         paste(reg$dataset, collapse = ", "), "\n(registry: ", WASP_DATASETS, ")")
  }
  cfg <- as.list(row)
  cfg$data_dir     <- file.path(REPO_DIR, "human", ds)
  cfg$analysis_dir <- file.path(cfg$data_dir, "wasp_refbias")
  cfg$fig_dir      <- file.path(DELIV_DIR, "figures",     ds)
  cfg$src_dir      <- file.path(DELIV_DIR, "source_data", ds)
  cfg$tab_dir      <- file.path(DELIV_DIR, "data",        ds)
  cfg$report_dir   <- file.path(DELIV_DIR, "report",      ds)
  cfg
}

# Set the module-level paths every figure section uses. Called once, at the top
# of the plotting script, from its dataset argument.
wasp_use_dataset <- function(ds) {
  cfg <- wasp_dataset(ds)
  DATASET     <<- ds
  DATASET_CFG <<- cfg
  LABEL       <<- cfg$label
  FIG_DIR     <<- cfg$fig_dir
  SRC_DIR     <<- cfg$src_dir
  DATA_DIR    <<- cfg$tab_dir
  REPORT_DIR  <<- cfg$report_dir

  # MINCOV is READ BACK from the tables, never hardcoded here. It used to be a
  # literal 10 in this file while WASP_prep_source_data.R built the evaluable
  # set at 20, so the bar panel's subtitle stated a threshold its own data did
  # not use. The generating stage records the value; the figure quotes it.
  MINCOV <<- wasp_mincov()

  for (d in c(FIG_DIR, SRC_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  message("dataset : ", ds, "  (", cfg$label, ")")
  message("mincov  : ", MINCOV, "  [read from report/", ds, "/WASP_summary.tsv]")
  invisible(cfg)
}

# Cross-dataset output. A figure that spans datasets cannot live in any one
# dataset's folder, so it gets its own. DATA_DIR / REPORT_DIR are deliberately
# NOT set: a cross-dataset script must read each dataset's tables explicitly
# via wasp_dataset(), so it can never silently plot one dataset's numbers
# under another's label.
wasp_use_combined <- function() {
  DATASET    <<- "combined"
  FIG_DIR    <<- file.path(DELIV_DIR, "figures",     "combined")
  SRC_DIR    <<- file.path(DELIV_DIR, "source_data", "combined")
  DATA_DIR   <<- NULL
  REPORT_DIR <<- NULL
  for (d in c(FIG_DIR, SRC_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  message("output  : figures/combined + source_data/combined")
  invisible(NULL)
}

# The coverage floor the CURRENT tables were built with.
wasp_mincov <- function() {
  f <- file.path(REPORT_DIR, "WASP_summary.tsv")
  if (!file.exists(f)) {
    stop("Missing ", f, "\nRun: Rscript ",
         file.path("scripts", "WASP_prep_source_data.R"), " ", DATASET)
  }
  s <- fread(f)
  v <- s[metric == "min_coverage", value]
  if (length(v) != 1L || is.na(suppressWarnings(as.integer(v)))) {
    stop("No usable min_coverage row in ", f)
  }
  as.integer(v)
}

# ── Figure writer ───────────────────────────────────────────────────────────
# Thin wrapper over save_panel(): SVG into figures/<dataset>/, every TSV into
# source_data/<dataset>/, so a figure folder never mixes panels with tables.
#
# `data = NULL` is for the combined figure, whose rows are exactly the rows
# already written by the panels it assembles; its sources= lists those
# per-panel source-data TSVs rather than duplicating them a second time.
save_figure <- function(plot, stem, width, height, sources, data = NULL,
                        grid = NA_character_, rasterized = FALSE,
                        script = NA_character_) {
  svg <- save_panel(plot, file.path(FIG_DIR, stem),
                    width = width, height = height,
                    sources = sources, data = data, grid = grid,
                    rasterized = rasterized, script = script,
                    sidecar_dir = SRC_DIR)

  message("  figure      : ", file.path("figures", DATASET, basename(svg)))
  if (is.null(data)) {
    message("  source data : ", length(sources), " panel tables (see ",
            file.path("source_data", DATASET, paste0(stem, ".sources.tsv")), ")")
  } else {
    message("  source data : ",
            file.path("source_data", DATASET, paste0(stem, ".source_data.tsv")),
            " (", format(nrow(as.data.table(data)), big.mark = ","), " rows)")
  }
  invisible(svg)
}
