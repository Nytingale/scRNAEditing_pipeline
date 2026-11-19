#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_hex_dedup_coord_vs_umi_plot.R   (2026-10-25; the hex section of GTDEF_bam_processing_plot.R)
# Pairs of de-duplication methods (coordinate, umi_tools, UMI-aware MarkDuplicates; 2026-11-16), one point per
# candidate site and tier (hex bins, rasterized):
#   GTDEF_hex_editlevel_coord_vs_umi        response figure 11: A-to-G editing level, per arm x tier
#   GTDEF_hex_depletion_ratio_coord_vs_umi  KO/NTC editing ratio, per KO arm x tier (ratios > 2 drawn at 2)
#   GTDEF_panel_hex_dedup_coord_vs_umi      the two, a over b
# Input : source_data/pm_gt/GTDEF_hex_dedup_coord_vs_umi.{el_points,dep_points}.tsv.gz, .{el_rho,dep_rho}.tsv
#         (GTDEF_hex_dedup_coord_vs_umi_source_data.R)
# Output: figures/pm_gt/GTDEF_{hex_editlevel,hex_depletion_ratio,panel_hex_dedup}_coord_vs_umi.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
SCRIPT <- "publication/gt_definition/scripts/GTDEF_hex_dedup_coord_vs_umi_plot.R"
STEM <- "GTDEF_hex_dedup_coord_vs_umi"

pts <- read_part(STEM, "el_points")
rho <- read_part(STEM, "el_rho")
pv <- pts[, .(p = spearman_p(edit_x, edit_y), n = .N), by = .(pair, tier, arm)]
rho <- merge(rho, pv, by = c("pair", "tier", "arm"), all.x = TRUE, sort = FALSE)
stopifnot(rho[!is.na(p), all(as.integer(gsub("\\D", "", sub(".*n = *", "", lab))) == n)])
rho[!is.na(p), lab := sub("^\u03c1 = ([^\n]*)\n", "Spearman's \u03c1 = \\1\n", lab)]
rho[!is.na(p), lab := unname(mapply(function(l, q) sub("\n", paste0("\n", q, "\n"), l, fixed = TRUE), lab, p_lab(p)))]
dep <- read_part(STEM, "dep_points")
rho_d <- read_part(STEM, "dep_rho")

hex_layer <- function(bins = 50, dpi = 900) {
  g <- geom_hex(bins = bins)
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::rasterize(g, dpi = dpi)
  } else {
    message("NOTE: ggrastr not installed; writing non-rasterized hex")
    g
  }
}

p_el <- ggplot(pts, aes(edit_x, edit_y)) +
  hex_layer(bins = 50) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.3) +
  facet_grid(arm ~ pair + tier) +
  geom_text(data = rho, aes(x = 0.02, y = 1.36, label = lab),
            hjust = 0, vjust = 1, size = 2.6, lineheight = 0.95, inherit.aes = FALSE) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1.36)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Editing level (x method)", y = "Editing level (y method)",
       title = "A-to-G editing at candidate sites by de-duplication method") +
  theme_manuscript(facets = TRUE) + grid_both +
  theme(panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(1.15, "lines"),
        legend.title.position = "top")
save_figure_sd(p_el, "GTDEF_hex_editlevel_coord_vs_umi", width = 17, height = 9, parts = c("el_points", "el_rho"),
               grid = "both", script = SCRIPT, rasterized = requireNamespace("ggrastr", quietly = TRUE),
               src_stem = STEM)

p_dep <- ggplot(dep, aes(ratio_x, ratio_y)) +
  hex_layer(bins = 60) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.3) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.3) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  facet_grid(ko_arm ~ pair + tier) +
  geom_text(data = rho_d, aes(x = 1.98, y = 0.08, label = lab),
            hjust = 1, vjust = 0, size = 3, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 2), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 2), oob = scales::squish) +
  coord_equal() +
  labs(x = "KO/NTC editing ratio (x method)", y = "KO/NTC editing ratio (y method)",
       title = "Depletion ratio by de-duplication method",
       caption = "Ratios above 2 are drawn at 2; dashed lines = 50% depletion") +
  theme_manuscript(facets = TRUE) + grid_both +
  theme(panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(1.15, "lines"),
        legend.title.position = "top")
save_figure_sd(p_dep, "GTDEF_hex_depletion_ratio_coord_vs_umi", width = 17, height = 5.5,
               parts = c("dep_points", "dep_rho"), grid = "both", script = SCRIPT,
               rasterized = requireNamespace("ggrastr", quietly = TRUE), src_stem = STEM)

p_panel <- cowplot::plot_grid(p_el + theme(plot.title = element_blank()),
                              p_dep + theme(plot.title = element_blank()),
                              ncol = 1, rel_heights = c(1.3, 1),
                              labels = "auto", label_size = 16, align = "v", axis = "lr")
save_figure_sd(p_panel, "GTDEF_panel_hex_dedup_coord_vs_umi", width = 17.5, height = 12.5,
               parts = c("el_points", "el_rho", "dep_points", "dep_rho"), script = SCRIPT, src_stem = STEM)
