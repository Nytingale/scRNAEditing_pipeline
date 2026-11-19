#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_config_factor_plot.R   (2026-10-15; draws from source data only since 2026-10-25)
# Response figure 13: paired F1 difference between two configurations that differ in one choice only, per
# tier and truth set (teratoma 6 KO + 3 REDIportal, held-out aorta 6 KO; v9, 76 configurations); one column per
# contrast (first column: each caller minus BCFtools), truth sets on the y-axis, points pool the tiers.
# Input : source_data/pm_gt/TERATOMA_v8_config_factor_delta_f1.pairs.tsv (TERATOMA_v8_config_factor_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_config_factor_delta_f1.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_config_factor_plot.R"
STEM <- "TERATOMA_v8_config_factor_delta_f1"

d <- read_part(STEM, "pairs")
x_max <- max(abs(d$delta)) * 1.05
# the SNP removal column label is one long line
# K scales every size (fonts, lines, points, spacing, page) together, so the figure keeps its look at any size
K <- 0.5
col_w <- ifelse(unique(d[order(Detail), .(Detail, Choice)])$Choice == "SNP removal", 1.4, 1)

p <- ggplot(d, aes(delta, Truth)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4 * K) +
  geom_point(aes(colour = method), position = position_jitter(width = 0, height = 0.22, seed = 1), size = 1.3 * K, alpha = 0.55) +
  geom_boxplot(fill = NA, colour = "black", width = 0.6, outlier.shape = NA, linewidth = 0.5 * K) +
  ggh4x::facet_nested(Family ~ Choice + Detail, scales = "free_y", space = "free_y",
                      strip = ggh4x::strip_nested(text_x = ggh4x::elem_list_text(face = c("bold", "plain")),
                                                  by_layer_x = TRUE)) +
  scale_colour_manual_cvd(values = cmap_method, name = NULL,
                          guide = guide_legend(nrow = 1, override.aes = list(size = 3 * K, alpha = 1))) +
  labs(x = "F1 difference (first minus second configuration)", y = NULL) +
  ggh4x::force_panelsizes(cols = col_w) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-x_max, x_max)) +
  theme_manuscript(base_size = 13 * K, facets = TRUE) + grid_both +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = -90), panel.spacing.x = unit(17 * K, "pt"),
        legend.key.size = unit(17 * K, "pt"), panel.grid.major = element_line(linewidth = 0.25 * K, colour = "grey85"))
save_figure_sd(p, STEM, width = 28 * K, height = 13 * K, parts = c("pairs", "summary"), grid = "both", script = SCRIPT)
