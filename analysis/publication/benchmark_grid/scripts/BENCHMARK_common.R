#!/usr/bin/env Rscript
# ============================================================================
# BENCHMARK_common.R
# Shared paths + figure-saving wrapper for the benchmark-grid deliverable.
# Sourced by BENCHMARK_prep_source_data.R and BENCHMARK_grid_source_data.R. Do not put figure code here.
# ============================================================================
# Naming convention for this deliverable (see ../../README.md):
#   BENCHMARK_<plottype>_<description>   one stem shared by the figure and its
#                                        source-data TSV.
# BENCHMARK_grid_plot.R builds every lettered panel and the combined figure they
# assemble into, from the tables of BENCHMARK_grid_source_data.R. The BENCHMARK_
# prefix groups them in a listing; figure numbers are assigned later from
# publication/FIGURE_INDEX.tsv.
#
# figures/ holds SVG only. save_figure() routes every TSV save_panel() emits
# (<stem>.sources.tsv, <stem>.source_data.tsv, figure_source_map.tsv) into
# source_data/ instead.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ── Locate the repo so this deliverable can be run from anywhere ────────────
bench_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) return(dirname(normalizePath(f[1], mustWork = FALSE)))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(o)) return(dirname(normalizePath(o, mustWork = FALSE)))
  getwd()
}

SCRIPT_DIR <- bench_script_dir()
DELIV_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_DIR   <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

source(file.path(REPO_DIR, "plotting", "plotting_functions.R"))

# ── EVERY PATH IS PER DATASET ───────────────────────────────────────────────
# This analysis runs on more than one dataset (per-sample, merged arms,
# Replogle, celltype-as-sample), so all four output trees are split by dataset
# with figure stems IDENTICAL across them — the folder says which dataset, so
# the same panel from two tiers is directly comparable and adding a dataset
# renames nothing. See publication/README.md.
#
# This was previously FLAT, which was a latent bug rather than a style issue:
# BENCHMARK_prep_source_data.R already took a DATASET argument, so a second
# dataset silently OVERWROTE the first's tables and figures.
#
# BENCHMARK_datasets.tsv is the single source of truth, shared with the bash
# stage, so a dataset is described once and the two cannot disagree.
BENCHMARK_DATASETS <- file.path(SCRIPT_DIR, "BENCHMARK_datasets.tsv")

bench_registry <- function() {
  if (!file.exists(BENCHMARK_DATASETS))
    stop("dataset registry not found: ", BENCHMARK_DATASETS)
  fread(BENCHMARK_DATASETS, colClasses = "character")
}

# Sets DATASET and the four per-dataset output dirs as globals, the same way
# WASP_common.R does, so the stage scripts stay simple.
bench_set_dataset <- function(ds) {
  reg <- bench_registry()
  if (!ds %in% reg$dataset)
    stop("unknown dataset '", ds, "'. Known: ", paste(reg$dataset, collapse = ", "),
         "\n(registry: ", BENCHMARK_DATASETS, ")")
  cfg <- reg[dataset == ds]

  DATASET     <<- ds
  DATASET_CFG <<- cfg
  DATASET_LABEL <<- cfg$label[1]
  # Resolved against the repo root so the registry can hold relative paths.
  SUMMARY_TSV <<- if (startsWith(cfg$summary_tsv[1], "/")) cfg$summary_tsv[1]
                  else file.path(REPO_DIR, cfg$summary_tsv[1])

  FIG_DIR    <<- file.path(DELIV_DIR, "figures",     ds)
  SRC_DIR    <<- file.path(DELIV_DIR, "source_data", ds)
  DATA_DIR   <<- file.path(DELIV_DIR, "data",        ds)
  REPORT_DIR <<- file.path(DELIV_DIR, "report",      ds)
  for (d in c(FIG_DIR, SRC_DIR, DATA_DIR, REPORT_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(cfg)
}

# ── The benchmark grid ──────────────────────────────────────────────────────
# The dimensions the pipeline varies. `sample` is the REPLICATE dimension: the
# panels here aggregate over it to a median with an IQR ribbon, and every other
# dimension stays distinguished by facet or aesthetic.
#
# >>> AGGREGATION OVER `sample` IS DELIBERATE AND IS DECLARED. <<<
# Collapsing a dimension on purpose and collapsing one by accident produce
# identical-looking output, so every panel passes aggregated = AGG_DIMS to
# assert_no_dropped_combos() and the guard reports what was collapsed and over
# how many levels. Anything neither shown nor declared warns. This is the
# standing rule for this project: in a grid of method x preprocessing x filter
# combinations, no unique combination may be silently dropped when collating or
# plotting. It exists because `snp_subtract` was silently absent from all six of
# the old in-pipeline panels, which fabricated line trajectories, hid bars
# behind bars and summed stacked bars to double height.
GRID_DIMS <- c("caller", "bam_type", "strand_mode", "snp_subtract")
AGG_DIMS  <- "sample"

# Filter-step order, shared by every panel so the x axis never disagrees.
STEP_LABELS <- c("1:multiMM", "2:noSNP", "3:stranded", "4:covEdit",
                 "5:noMM(w1)", "6:clustered(w2)", "7:no100pct")

# ── Figure writer ───────────────────────────────────────────────────────────
# Thin wrapper over save_panel(): SVG into figures/, every TSV into
# source_data/, so figures/ never mixes rendered panels with tables.
#
# `data = NULL` is for a combined figure whose rows are exactly the rows already
# written by the panels it assembles; its sources= lists those per-panel
# source-data TSVs rather than duplicating them.
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
