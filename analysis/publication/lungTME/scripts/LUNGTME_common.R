#!/usr/bin/env Rscript
# ============================================================================
# LUNGTME_common.R
# Shared paths + figure-saving wrapper for the lung tumour microenvironment application
# (Reviewer 1, Major 5; Reviewer 3, Major 10). Analysis outputs live under human/dl/;
# this deliverable holds only figures and their source data. Do not put figure code here.
#
# Naming convention (see ../../README.md):
#   LUNGTME_<plottype>_<description>
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

lungtme_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) return(dirname(normalizePath(f[1], mustWork = FALSE)))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(o)) return(dirname(normalizePath(o, mustWork = FALSE)))
  getwd()
}

SCRIPT_DIR <- lungtme_script_dir()
DELIV_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_DIR   <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

source(file.path(REPO_DIR, "plotting", "plotting_functions.R"))

lungtme_set_dataset <- function(ds = "dl") {
  DATASET <<- ds
  FIG_DIR <<- file.path(DELIV_DIR, "figures",     ds)
  SRC_DIR <<- file.path(DELIV_DIR, "source_data", ds)
  for (d in c(FIG_DIR, SRC_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(ds)
}

DEPMAP_RES <- file.path(REPO_DIR, "human", "dl", "depmap_results")

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
