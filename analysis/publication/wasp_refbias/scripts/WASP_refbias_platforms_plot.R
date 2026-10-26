#!/usr/bin/env Rscript
# ============================================================================
# WASP_refbias_platforms_plot.R   (split from WASP_refbias_plot_platforms.R 2026-10-25; source data only)
# WASP reference bias, both platforms and all coverage floors; response figure 15 = WASP_panel_refbias_platforms:
#   a  WASP_hist_delta_editing_level_platforms  the correction is negligible (per-site editing level shift)
#   b  WASP_bar_removed_read_alleles_platforms  it costs ~88% of reads at the sites (gap = bias, CI clustered by site)
#   c  WASP_bar_site_attrition_platforms        so ~90% of measurable sites go
# (a) is one row of 6 facets (platform x floor); (b) and (c) are grouped bars, platform on x, faceted by floor, so
# the facet columns do not align between panels; align = "v" keeps the plot areas and y axes registered.
# Input : source_data/combined/WASP_refbias_platforms.<part>.tsv (WASP_refbias_platforms_source_data.R)
# Output: figures/combined/WASP_{hist_delta_editing_level,bar_removed_read_alleles,bar_site_attrition,
#         panel_refbias}_platforms.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "WASP_plot_common.R"))
SCRIPT <- "publication/wasp_refbias/scripts/WASP_refbias_platforms_plot.R"
STEM <- "WASP_refbias_platforms"

dat.a <- read_part(STEM, "a_delta")
stats.a <- read_part(STEM, "a_stats")
dat.b <- read_part(STEM, "b_removed")
gaps.b <- read_part(STEM, "b_gaps")
dat.c <- read_part(STEM, "c_sites")
drops.c <- read_part(STEM, "c_drops")

# The binding constraint is panel (a)'s 3-line stats block, which has only ~2.5 in of facet width; 4.0 still clears it.
ANNOT_SIZE <- 4.0
BASE_SIZE <- 13
TOP_PAD_LINES <- 0.75
TOP_VJUST_2LINE <- 1 + TOP_PAD_LINES / 2
TOP_VJUST_3LINE <- 1 + TOP_PAD_LINES / 3
TOP_EXPAND <- expansion(mult = c(0, 0.42))
DODGE <- position_dodge(width = 0.72)

p.a <- ggplot(dat.a, aes(dEL)) +
  geom_histogram(bins = 50, fill = unname(okabe_ito["blue"]), colour = "white", linewidth = 0.1) +
  geom_text(data = stats.a, aes(x = -1, y = Inf, label = lab),
            hjust = 0, vjust = TOP_VJUST_3LINE, size = ANNOT_SIZE, lineheight = 1.05, inherit.aes = FALSE) +
  # facet_nested with a free y per facet: the evaluable set shrinks ~4x across floors and differs between platforms
  ggh4x::facet_nested(~ platform + cov, scales = "free_y") +
  scale_y_continuous(labels = shorten_axis, expand = TOP_EXPAND) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(title = "Per-site editing level shifts",
       subtitle = "The correction WASP applies is negligible on both platforms, at every coverage floor",
       x = expression(Delta * " editing level (WASP-passed - all reads)"), y = "Number of sites") +
  theme_manuscript(base_size = BASE_SIZE, facets = TRUE) + grid_y +
  theme(panel.spacing.x = unit(1.2, "lines"))
save_figure_sd(p.a, "WASP_hist_delta_editing_level_platforms", width = 15.0, height = 5.0,
               parts = c("a_delta", "a_stats"), grid = "y", script = SCRIPT, src_stem = STEM)

p.b <- ggplot(dat.b, aes(platform, removal_rate, fill = allele)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.2, position = DODGE) +
  geom_text(aes(label = sprintf("%.2f%%", 100 * removal_rate)), position = DODGE, vjust = -0.35, size = ANNOT_SIZE) +
  geom_text(data = gaps.b, aes(x = 1.5, y = Inf, label = lab),
            vjust = TOP_VJUST_2LINE, size = ANNOT_SIZE, lineheight = 1.15, inherit.aes = FALSE) +
  facet_wrap(~ cov, nrow = 1) +
  scale_fill_manual_cvd(values = cmap_allele, name = NULL) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = TOP_EXPAND) +
  # neutral subtitle: the clustered CI spans zero for jzlxv3 and not for jzlss2; each facet states its own
  labs(title = "Reads removed at REDIportal sites",
       subtitle = "Bar height = the cost; the gap within each pair = the bias (95% CI clustered by site)",
       x = NULL, y = "Percentage of allele's reads removed") +
  theme_manuscript(base_size = BASE_SIZE, facets = TRUE) + grid_y +
  theme(panel.spacing.x = unit(1.2, "lines"))
save_figure_sd(p.b, "WASP_bar_removed_read_alleles_platforms", width = 15.0, height = 5.0,
               parts = c("b_removed", "b_gaps"), grid = "y", script = SCRIPT, src_stem = STEM)

# shared y across facets: the covered set itself shrinks as the floor rises, part of the cost
p.c <- ggplot(dat.c, aes(platform, n_sites, fill = set)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.2, position = DODGE) +
  geom_text(aes(label = unlist(lapply(dat.c$n_sites, function(x) shorten_num(x)))),
            position = DODGE, vjust = -0.35, size = ANNOT_SIZE) +
  geom_text(data = drops.c, aes(x = 1.5, y = Inf, label = lab),
            vjust = TOP_VJUST_2LINE, size = ANNOT_SIZE, lineheight = 1.15, inherit.aes = FALSE) +
  facet_wrap(~ cov, nrow = 1) +
  scale_fill_manual_cvd(values = cmap_attrition, name = NULL) +
  scale_y_continuous(labels = shorten_axis, expand = TOP_EXPAND) +
  labs(title = "REDIportal sites measurable before and after WASP",
       subtitle = "~90% of the sites where editing can be measured are lost to buy that correction",
       x = NULL, y = "Number of A-to-I sites") +
  theme_manuscript(base_size = BASE_SIZE, facets = TRUE) + grid_y +
  theme(panel.spacing.x = unit(1.2, "lines"))
save_figure_sd(p.c, "WASP_bar_site_attrition_platforms", width = 15.0, height = 5.0,
               parts = c("c_sites", "c_drops"), grid = "y", script = SCRIPT, src_stem = STEM)

p.panel <- cowplot::plot_grid(p.a, p.b, p.c, ncol = 1, align = "v", axis = "lr", labels = "auto", label_size = 20)
save_figure_sd(p.panel, "WASP_panel_refbias_platforms", width = 12.0, height = 15.0,
               parts = c("a_delta", "a_stats", "b_removed", "b_gaps", "c_sites", "c_drops"), grid = "y",
               script = SCRIPT, src_stem = STEM)
