#!/usr/bin/env Rscript
# ============================================================================
# LUNGTME_depmap_exADAR_plot.R   (2026-10-04; draws from source data only since 2026-10-25)
# Response figure 19 (Reviewer 3, Major 10). DepMap 26Q1 lung cell lines; ISG scores exclude only ADAR.
# Dependent = DepMap gene dependency probability >= 0.5.
#   a  ADAR and ADARB1 dependency probability, all lung lines; dependent lines in red
#   b  per-gene Spearman's rho (expression vs ADAR Chronos; one point = one gene), dependent vs
#      non-dependent lines, one facet per ISG set plus all expressed genes as background (rasterised);
#      p = paired Wilcoxon signed-rank test, dependent vs non-dependent rho of the same genes
#   c  single genes vs ADAR Chronos (Spearman): the dsRNA sensors and the 20 ISGs with the most
#      negative rho in dependent lines; non-dependent lines shown for comparison; dashed line = cut-off
#      of the 5% lowest rho values over all expressed genes in dependent lines
# Input : source_data/dl/LUNGTME_panel_depmap_adar_isg_exADAR.<part>.tsv (LUNGTME_depmap_exADAR_source_data.R)
# Output: figures/dl/LUNGTME_panel_depmap_adar_isg_exADAR.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "LUNGTME_plot_common.R"))
SCRIPT <- "publication/lungTME/scripts/LUNGTME_depmap_exADAR_plot.R"
STEM <- "LUNGTME_panel_depmap_adar_isg_exADAR"
PARTS <- c("a_prob", "b_rho", "b_tests", "c_genes", "c_ref")

la <- read_part(STEM, "a_prob")
lb <- read_part(STEM, "b_rho")
p_b_test <- read_part(STEM, "b_tests")
lc <- read_part(STEM, "c_genes")
lab_q <- read_part(STEM, "c_ref")

PROB_CUT <- 0.5
cmap_dep <- c(Dependent = "#CB181D", `Not dependent` = "grey55")
LAB_SIZE <- 5
jit <- position_jitter(width = 0.2, height = 0, seed = 1)
base <- theme_manuscript(base_size = 18, facets = TRUE)

# a
p_a <- ggplot(la, aes(gene, prob)) +
  geom_hline(yintercept = PROB_CUT, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = status), position = jit, size = 1.6, alpha = 0.7) +
  geom_boxplot(fill = NA, outlier.shape = NA, colour = "black", width = 0.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.03, 0.03))) +
  scale_colour_manual(values = cmap_dep, name = NULL) +
  labs(x = NULL, y = "Dependency probability") + base + grid_y + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# b
y_top <- lb[, max(rho)] + 0.05
p_b <- ggplot(lb, aes(status, rho)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
  geom_point(data = lb[bg == FALSE], aes(colour = status), position = jit, size = 0.9, alpha = 0.45) +
  ggrastr::rasterise(geom_point(data = lb[bg == TRUE], aes(colour = status), position = jit, size = 0.3, alpha = 0.12),
                     dpi = 600) +
  geom_boxplot(fill = NA, outlier.shape = NA, colour = "black", width = 0.55) +
  geom_segment(data = p_b_test, aes(x = 1, xend = 2, y = y_top, yend = y_top), inherit.aes = FALSE) +
  geom_text(data = p_b_test, aes(x = 1.5, y = y_top, label = lab), vjust = -0.5, size = LAB_SIZE,
            inherit.aes = FALSE) +
  scale_colour_manual(values = cmap_dep, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.15))) +
  facet_wrap(~strip, nrow = 1) +
  labs(x = NULL, y = "Spearman's rho, expression vs\ngrowth effect of ADAR KO") +
  base + grid_y + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text = element_text(size = 14))

# c
REF_COL <- "#2171B5"
p_c <- ggplot(lc, aes(gene, rho)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
  geom_hline(yintercept = lab_q$rho, colour = REF_COL, linetype = "dashed") +
  geom_text(data = lab_q, aes(label = lab), hjust = 1, vjust = -0.5, size = LAB_SIZE + 0.5, colour = REF_COL) +
  geom_point(aes(colour = status), size = 2.6) +
  scale_colour_manual(values = cmap_dep, name = NULL) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Spearman's rho, expression vs\ngrowth effect of ADAR KO") +
  base + grid_y + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), legend.position = "right")

p <- patchwork::wrap_plots(patchwork::wrap_plots(p_a, p_b, widths = c(0.2, 1)),
                           p_c, ncol = 1, heights = c(1, 0.9)) +
  patchwork::plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 18, face = "bold"))
save_figure_sd(p, STEM, width = 16, height = 11.5, parts = PARTS, rasterized = TRUE, script = SCRIPT)
