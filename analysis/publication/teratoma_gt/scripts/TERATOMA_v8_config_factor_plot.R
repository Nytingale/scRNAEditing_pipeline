#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_config_factor_plot.R   (2026-10-15; draws from source data only since 2026-10-25)
# Response figure 10: paired F1 difference between two configurations that differ in one choice only, per
# tier and v8 truth set (6 KO + 3 REDIportal); one column per contrast, truth sets on the y-axis, points pool the
# three tiers.
# Input : source_data/pm_gt/TERATOMA_v8_config_factor_delta_f1.pairs.tsv (TERATOMA_v8_config_factor_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_config_factor_delta_f1.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_config_factor_plot.R"
STEM <- "TERATOMA_v8_config_factor_delta_f1"

d <- read_part(STEM, "pairs")

p <- ggplot(d, aes(delta, Truth)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(aes(colour = method), position = position_jitter(width = 0, height = 0.22, seed = 1), size = 1.3, alpha = 0.55) +
  geom_boxplot(fill = NA, colour = "black", width = 0.6, outlier.shape = NA, linewidth = 0.5) +
  ggh4x::facet_nested(Family ~ Choice + Detail, scales = "free_y", space = "free_y",
                      strip = ggh4x::strip_nested(text_x = ggh4x::elem_list_text(face = c("bold", "plain")),
                                                  by_layer_x = TRUE)) +
  scale_colour_manual_cvd(values = cmap_method, name = NULL,
                          guide = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 1))) +
  labs(x = "F1 difference (first minus second configuration)", y = NULL,
       title = "Pairs of configuration sets differing in one choice") +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  theme_manuscript(base_size = 13, facets = TRUE) + grid_both +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"))
save_figure_sd(p, STEM, width = 19, height = 7.5, parts = c("pairs", "summary"), grid = "both", script = SCRIPT)
