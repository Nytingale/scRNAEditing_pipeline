# ============================================================================
# AORTA_plot_common.R   (2026-10-27)
# Sourced by the AORTA_*_source_data.R and AORTA_*_plot.R scripts. The plot scripts draw a figure from its source
# data only (source_data/aorta_gt/<stem>.<part>.tsv, written by the matching *_source_data.R script).
# Aorta: GEO GSE255661, Myh11+ cells, 3 Adar1 conditional-KO vs 3 control mice (mouse/aorta_gse255661/).
# Teratoma: v8 truth (human/pm/gt_eval_v8/). Replogle is not part of the aorta comparison.
# ============================================================================
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
DELIV_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_DIR   <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
source(file.path(REPO_DIR, "plotting", "plotting_functions.R"))
source(file.path(REPO_DIR, "plotting", "gt_v8_paths.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))

DATASET <- "aorta_gt"
FIG_DIR <- file.path(DELIV_DIR, "figures", DATASET)
SRC_DIR <- file.path(DELIV_DIR, "source_data", DATASET)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
AO_DIR  <- file.path(REPO_DIR, "mouse", "aorta_gse255661")
AO_EVAL <- file.path(AO_DIR, "gt_eval_v8")

COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
DATASET_LAB <- c(teratoma = "Teratoma (human)", aorta = "Aorta (mouse)")

read_part <- function(stem, part) read_sd(SRC_DIR, stem, part)
save_figure_sd <- function(plot, stem, width, height, parts, grid = NA_character_, script = NA_character_,
                           src_stem = stem) {
  save_panel(plot, file.path(FIG_DIR, stem), width = width, height = height,
             sources = sd_sources(SRC_DIR, src_stem, parts), grid = grid, script = script, sidecar_dir = SRC_DIR)
}
cfg_label <- function(d) sprintf("%s | %s | %s | %s", relabel(d$caller, label_caller), relabel(d$bam_type, label_bamType),
                                 relabel(d$strand_mode, label_strandMode), relabel(d$snp_subtract, label_snpSubtract))
