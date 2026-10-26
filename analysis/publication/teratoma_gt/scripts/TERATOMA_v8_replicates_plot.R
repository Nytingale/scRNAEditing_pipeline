#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_replicates_plot.R   (2026-10-15; draws from source data only since 2026-10-25)
# Replicate-aware calling vs pooled calling (review R3): JACUSA2 replicate-aware vs pooled F1, one point per
# configuration x truth; the 4 JACUSA2 configurations on deduplicated BAMs (shape = preprocessing, as
# TERATOMA_v8_overall_vs_ag_plot.R); recall over all positives (pooled-arm coverage) and over positives with
# >= 10 reads in each of T1-T4.
# Input : source_data/pm_gt/TERATOMA_v8_replicates.f1_pairs.tsv (TERATOMA_v8_replicates_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_replicates.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_replicates_plot.R"
STEM <- "TERATOMA_v8_replicates"

b <- read_part(STEM, "f1_pairs")
FAM_LAB <- levels(b$fam)
PROC_LEVELS <- levels(b$proc)
shape_proc <- setNames(c(21, 24), PROC_LEVELS)

p_b <- ggplot(b, aes(pooled, jacusa2rep, fill = fam, shape = proc, size = proc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(colour = "grey15", stroke = 0.3, alpha = 0.7) +
  scale_size_manual(values = setNames(c(4.2, 2.6), PROC_LEVELS), guide = "none") +
  facet_wrap(~ Set, nrow = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_fill_manual_cvd(values = setNames(c("#3182BD", "#E6550D"), FAM_LAB), name = NULL,
                        guide = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, size = 4, alpha = 0.9))) +
  scale_shape_manual(values = shape_proc, name = NULL,
                     guide = guide_legend(nrow = 1, order = 2, override.aes = list(fill = "grey60", size = 4, alpha = 0.9))) +
  labs(x = "F1, pooled JACUSA2", y = "F1, replicate-aware JACUSA2", title = NULL,
       subtitle = "One point per JACUSA2 configuration x ground truth; deduplicated BAMs") +
  theme_manuscript(base_size = 12, facets = TRUE) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(2, "pt"))

save_figure_sd(p_b, STEM, width = 10, height = 6.5, parts = "f1_pairs", script = SCRIPT)
