#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_caller_resources_plot.R   (2026-10-15; draws from source data only since 2026-10-25)
# Response figure 12: run time, peak memory and storage of each variant calling task, per caller and tier.
# Input : source_data/pm_gt/TERATOMA_v8_caller_resources.tasks.tsv (TERATOMA_v8_caller_resources_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_caller_resources.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_caller_resources_plot.R"
STEM <- "TERATOMA_v8_caller_resources"

long <- read_part(STEM, "tasks")
p <- ggplot(long, aes(value, method, colour = method)) +
  geom_point(position = position_jitter(width = 0, height = 0.15, seed = 1), size = 0.8, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.5, colour = "black", fill = NA) +
  facet_grid(Tier ~ Metric, scales = "free_x") +
  ggh4x::facetted_pos_scales(x = list(
    Metric == "Run time (min)" ~ scale_x_log10(labels = function(x) format(x, drop0trailing = TRUE, scientific = FALSE)),
    Metric == "Peak memory (GB)" ~ scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))),
    Metric == "Storage (MB)" ~ scale_x_log10(labels = function(x) format(x, drop0trailing = TRUE, scientific = FALSE)))) +
  scale_colour_manual_cvd(values = cmap_method, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_manuscript(base_size = 13, facets = TRUE) + grid_x + theme(strip.text.y = element_text(angle = 0), panel.spacing.x = unit(0.7, "cm"))
save_figure_sd(p, STEM, width = 12.5, height = 8, parts = c("tasks", "summary"), grid = "x", script = SCRIPT)
