#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_replicates_arm_rank_plot.R   (2026-10-15; draws from source data only since 2026-10-25)
# Response figure 17: replicate-aware vs pooled JACUSA2 at the dataset (arm) tier.
#   a  replicate-aware vs pooled JACUSA2 F1, the 2 configurations on UMI-deduplicated BAMs x 9 truths
#   b  dataset-tier F1 of the 58 ranked configurations (small points, caller colour) and of replicate-aware JACUSA2
#      (large diamonds, UMI-deduplicated, jittered among them)
# Input : source_data/pm_gt/TERATOMA_v8_replicates_arm_rank.{a_f1_pairs,b_pooled,b_replicate}.tsv
#         (TERATOMA_v8_replicates_arm_rank_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_replicates_arm_rank.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_replicates_arm_rank_plot.R"
STEM <- "TERATOMA_v8_replicates_arm_rank"

sc <- read_part(STEM, "a_f1_pairs")
pool <- read_part(STEM, "b_pooled")
rep <- read_part(STEM, "b_replicate")

p_a <- ggplot(sc, aes(pooled, jacusa2rep)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(shape = 21, size = 3.6, fill = cmap_method[["JACUSA2"]], colour = "grey15", stroke = 0.3, alpha = 0.8) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "F1, pooled JACUSA2", y = "F1, replicate-aware JACUSA2") +
  theme_manuscript(base_size = 12) + grid_both

REP_LAB <- "Replicate-aware JACUSA2\n(de-duplicated (UMI))"
p_b <- ggplot(mapping = aes(f1, TruthY)) +
  geom_point(data = pool, aes(colour = method),
             position = position_jitter(width = 0, height = 0.16, seed = 1), size = 1.3, alpha = 0.6) +
  geom_point(data = rep[, Rep := REP_LAB], aes(shape = Rep), fill = cmap_method[["JACUSA2"]], colour = "black",
             position = position_jitter(width = 0, height = 0.16, seed = 2), stroke = 0.7, size = 3) +
  scale_colour_manual_cvd(values = cmap_method, name = "Ranked configurations",
                          guide = guide_legend(nrow = 2, order = 1, title.position = "top",
                                               override.aes = list(size = 3, alpha = 1))) +
  scale_shape_manual(values = setNames(23, REP_LAB), name = NULL, guide = guide_legend(order = 3)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "F1, dataset tier", y = NULL) +
  theme_manuscript(base_size = 12) + grid_x

p <- (p_a | p_b) + plot_layout(widths = c(1, 1.6), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"), legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(face = "bold"), legend.justification = "top")
save_figure_sd(p, STEM, width = 13, height = 6.5, parts = c("a_f1_pairs", "b_pooled", "b_replicate"), script = SCRIPT)
