#!/usr/bin/env Rscript
# =============================================================================
# TERATOMA_cross_tier_grid_plot.R
#
# Tile grid of every reported configuration x mismatch-window-filter state, ranked
# by cross-tier F1 rank-sum (50% depletion). Local mismatch context (HANDOFF 2026-09-29
# UPDATE-22, UPDATE-24): near other MM types excluded, i.e. sites with other mismatch
# types nearby removed (step 5); + clustered, only sites with other candidate sites nearby
# kept (step 6, the benchmark call set). No-filter rows (step 4) are not shown
# (GRID_WINDOW_STATES). Ranks are over all rows of a tier.
# Version with the None rows: ../_superseded/drop_none_2026-09-30/scripts/
# Columns: Preprocessing | Caller | Stranding | SNP removal | Local mismatch context | Rank.
# Pre-window-column version: ../_superseded/window_column_2026-09-30/scripts/
#
# Usage: Rscript TERATOMA_cross_tier_grid_plot.R [DATASET]
# Alternative truth (HANDOFF UPDATE-25): GRID_SCORING_TSV=<gt_tier_scoring-format table>
#   GRID_STEM=<figure stem> GRID_TITLE_SUFFIX=<text> draws the same grid from that table.
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "TERATOMA_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
teratoma_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

# Data (50% depletion) ----
ALT_TSV <- Sys.getenv("GRID_SCORING_TSV")
SOURCES <- if (nzchar(ALT_TSV)) ALT_TSV else file.path(WINDOW_SCORING_DIRS[GRID_WINDOW_STATES], "gt_tier_scoring.tsv")
sc <- if (nzchar(ALT_TSV)) fread(ALT_TSV)[, window_filter := "both"] else
  read_window_scoring("gt_tier_scoring.tsv", GRID_WINDOW_STATES)
sc[, f1 := fifelse(precision + recall_vs_union > 0,
                   2 * precision * recall_vs_union / (precision + recall_vs_union),
                   NA_real_)]
sc[, f1_rank := frank(-f1, ties.method = "min"), by = tier]

rank_sum <- sc[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE),
                   method = relabel(caller[1], label_caller)),
               by = .(caller, bam_type, strand_mode, snp_subtract, window_filter)]
# ties on the rank sum are broken by the higher mean F1 across tiers
rank_sum[, overall_rank := frank(list(rank_sum, -mean_f1), ties.method = "min")]
setorder(rank_sum, overall_rank)

# Manuscript labels
rank_sum[, Preprocessing    := relabel(bam_type, label_bamType)]
rank_sum[, Caller           := relabel(caller, label_caller)]
rank_sum[, Stranding        := relabel(strand_mode, label_strandMode)]
rank_sum[, `SNP removal`    := relabel(snp_subtract, label_snpSubtract)]
rank_sum[, `Local mismatch context` := relabel(window_filter, label_windowFilter)]

rank_sum[, row_idx := .I]
rank_sum[, row := factor(row_idx, levels = rev(seq_len(.N)))]

# Colour palettes sourced from plotting_functions.R: cmap_proc, cmap_strand,
# cmap_snp, cmap_method, cmap_window.

# 5 main columns ----
cols <- c("Preprocessing", "Caller", "Stranding", "SNP removal", if (GRID_SHOW_WINDOW) "Local mismatch context")
id_cols <- c("row", "row_idx", "overall_rank", "rank_sum")
long <- melt(rank_sum[, c(id_cols, cols), with = FALSE],
             id.vars = id_cols, measure.vars = cols,
             variable.name = "column", value.name = "value")
long[, column := factor(column, levels = cols)]

# Assign fill
long[, fill := NA_character_]
long[column == "Preprocessing",  fill := cmap_proc[value]]
long[column == "Caller",         fill := cmap_method[value]]
long[column == "Stranding",      fill := cmap_strand[value]]
long[column == "SNP removal",    fill := cmap_snp[value]]
long[column == "Local mismatch context", fill := cmap_window[value]]

# Text colour: white on dark backgrounds, black on light.
# Threshold 120 so BCFtools blue (#00AAFF, lum 129) gets black text while
# GATK red (63), JACUSA2 purple (87), scAllele ochre (91) get white.
text_col <- function(hex) {
  rgb <- col2rgb(hex)
  lum <- rgb[1, ] * 0.299 + rgb[2, ] * 0.587 + rgb[3, ] * 0.114
  ifelse(lum < 120, "white", "black")
}
long[, txt_col := text_col(fill)]

# Main tile grid (5 columns) ----
p_main <- ggplot(long, aes(column, row)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.3) +
  scale_fill_identity() +
  geom_text(aes(label = value, colour = txt_col), size = 3.6, show.legend = FALSE) +
  scale_colour_identity() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(5, 0, 5, 5))

# Narrow rank strip ----
rank_dt <- rank_sum[, .(row, rank_label = as.character(overall_rank))]

p_rank <- ggplot(rank_dt, aes(x = factor(1), y = row)) +
  geom_tile(fill = "white", colour = "white", linewidth = 0.3) +
  geom_text(aes(label = rank_label), size = 3.4) +
  scale_x_discrete(position = "top", labels = "Rank") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(5, 5, 5, 0))

# Compose ----
title_grob <- cowplot::ggdraw() +
  cowplot::draw_label(paste0("All ", nrow(rank_sum),
                             " configurations", if (GRID_SHOW_WINDOW) " x local mismatch context",
                             " ranked by sum of cross-tier F1 ranks", Sys.getenv("GRID_TITLE_SUFFIX")),
                      fontface = "bold", size = 13, x = 0.5, hjust = 0.5)
grid_row <- cowplot::plot_grid(p_main, p_rank, ncol = 2,
                               rel_widths = c(5, 0.4), align = "h", axis = "tb")
fig_h <- 0.8 + 0.30 * nrow(rank_sum)   # one-line labels
p_final <- cowplot::plot_grid(title_grob, grid_row, ncol = 1,
                              rel_heights = c(0.5, fig_h - 0.5))

STEM <- Sys.getenv("GRID_STEM", "TERATOMA_cross_tier_grid")
save_figure(p_final, STEM, width = 10, height = fig_h,
            sources = SOURCES, data = rank_sum, grid = "none",
            script = SCRIPT)

message("Done -> figures/", DATASET, "/", STEM, ".svg")
