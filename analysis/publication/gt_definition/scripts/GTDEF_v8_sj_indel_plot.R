#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_sj_indel_plot.R   (2026-10-23; panels redone 2026-10-24; draws from source data only since 2026-10-25)
# Response figure 13 (R2 major 3): KO ground truth sites (>= 50% reduction) by distance, 0-10 nt and > 10 nt
# (own y scale), counts above the bars.
#   a  distance to the nearest splice junction
#   b  distance to the nearest indel in the merged NTC reads
# Input : source_data/pm_gt/GTDEF_v8_sj_indel.distances.tsv (GTDEF_v8_sj_indel_source_data.R)
# Output: figures/pm_gt/GTDEF_v8_sj_indel.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/gt_definition/scripts/GTDEF_v8_sj_indel_plot.R"
STEM <- "GTDEF_v8_sj_indel"

n <- read_part(STEM, "distances")
BS <- 14
FEAT_LAB <- c(junction = "Distance to nearest splice junction (nt)", indel = "Distance to nearest indel (nt)")
dodge <- position_dodge(width = 0.8)

panel <- function(feat) {
  ggplot(n[feature == feat], aes(d, n_sites, fill = k)) +
    geom_col(position = dodge, width = 0.75) +
    geom_text(aes(label = scales::label_comma()(n_sites), colour = k), position = dodge, vjust = -0.35,
              size = 4.2, show.legend = FALSE) +
    ggh4x::facet_grid2(g ~ rng, scales = "free", independent = "y", switch = "y") +
    ggh4x::force_panelsizes(cols = c(11, 2.2)) +
    scale_fill_manual_cvd(values = cmap_truth_class, name = "KO ground truth sites (\u2265 50% reduction)") +
    scale_colour_manual_cvd(values = cmap_truth_class, guide = "none") +
    scale_y_continuous(labels = shorten_axis, expand = expansion(mult = c(0, 0.16))) +
    labs(x = FEAT_LAB[[feat]], y = "Ground truth sites") +
    theme_manuscript(base_size = BS, facets = TRUE) + strip_y_outside + grid_y +
    theme(strip.text.x = element_blank(), legend.position = "top") +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
}

p <- (panel("junction") / panel("indel")) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top", plot.tag = element_text(face = "bold", size = BS + 2))
save_figure_sd(p, STEM, width = 8, height = 11, parts = "distances", script = SCRIPT)
