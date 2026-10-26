# ============================================================================
# WASP_plot_common.R
# Sourced by the WASP_*_plot.R scripts, which draw a figure from its source data only
# (source_data/<dataset>/<stem>.<part>.tsv, written by the matching *_source_data.R script).
# ============================================================================
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
DELIV_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_DIR   <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
source(file.path(REPO_DIR, "plotting", "plotting_functions.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))

# The cross-platform figures live in figures/combined and source_data/combined; the per-dataset ones in
# figures/<dataset> and source_data/<dataset> (wasp_plot_dataset()).
wasp_plot_dataset <- function(ds) {
  DATASET <<- ds
  FIG_DIR <<- file.path(DELIV_DIR, "figures", ds)
  SRC_DIR <<- file.path(DELIV_DIR, "source_data", ds)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
}
wasp_plot_dataset("combined")

read_part <- function(stem, part) read_sd(SRC_DIR, stem, part)

# src_stem: the source data stem when several figures are drawn from one set of parts.
save_figure_sd <- function(plot, stem, width, height, parts, grid = NA_character_, rasterized = FALSE,
                           script = NA_character_, src_stem = stem) {
  save_panel(plot, file.path(FIG_DIR, stem), width = width, height = height,
             sources = sd_sources(SRC_DIR, src_stem, parts), grid = grid, rasterized = rasterized,
             script = script, sidecar_dir = SRC_DIR)
}
