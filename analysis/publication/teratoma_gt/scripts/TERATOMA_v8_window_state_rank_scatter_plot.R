#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_window_state_rank_scatter_plot.R   (2026-10-18; draws from source data only since 2026-10-25)
# Response figure 7: overall cross-tier rank of the 58 ranked configurations without the discordant-mismatch filter
# against the main analysis, one panel per t = 50% KO truth. Axes reversed so rank 1 is top right.
# Input : source_data/pm_gt/TERATOMA_v8_window_state_rank_scatter.{ranks,top_k,rho}.tsv
#         (TERATOMA_v8_window_state_rank_scatter_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_window_state_rank_scatter.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_window_state_rank_scatter_plot.R"
STEM <- "TERATOMA_v8_window_state_rank_scatter"

w <- read_part(STEM, "ranks")
top_k <- read_part(STEM, "top_k")
st <- read_part(STEM, "rho")
pv <- w[, .(rho2 = cor(none, mismatch, method = "spearman"), p = spearman_p(none, mismatch)), by = Truth]
st <- merge(st, pv, by = "Truth")
stopifnot(st[, all(abs(rho - rho2) < 1e-9)])
shape_proc <- setNames(c(21, 24, 22, 23, 25)[seq_along(levels(w$proc))], levels(w$proc))
# 76 configurations: the shapes of TRUTHFREE_v10_harmonic_scatter.R
if ("De-duplicated (UMI-aware MarkDuplicates)" %in% levels(w$proc))
  shape_proc <- setNames(c(21, 22, 23, 24, 25, 4), c("De-duplicated (coordinate)", "De-duplicated (umi_tools)",
    "De-duplicated (UMI-aware MarkDuplicates)", "Recalibrated", "Recalibrated (Q20)", "Recalibrated (default)"))[levels(w$proc)]

N_CFG <- w[, uniqueN(config)]
BR <- unique(c(1, seq(10, N_CFG - 5, 10), N_CFG))
assert_cvd_safe(cmap_method, "scale_fill_manual")
p <- ggplot(w, aes(none, mismatch)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_rect(data = top_k, aes(xmin = -0.6, xmax = k + 1.6, ymin = -0.6, ymax = k + 1.6), inherit.aes = FALSE,
            fill = NA, colour = "black", linewidth = 0.5) +
  geom_point(aes(fill = method, shape = proc), colour = "grey15", stroke = 0.3, size = 3, alpha = 0.8) +
  geom_text(data = st, aes(x = 1, y = N_CFG, label = paste0(spearman_lab(rho), "\n", p_lab(p))), hjust = 1, vjust = 0, size = 4.4, lineheight = 0.95) +
  facet_wrap(~ Truth, nrow = 1) +
  scale_fill_manual(values = cmap_method, name = NULL,
                    guide = guide_legend(nrow = 2, order = 1, override.aes = list(shape = 21, size = 4.5))) +
  scale_shape_manual(values = shape_proc, name = NULL,
                     guide = guide_legend(nrow = 2, order = 2, override.aes = list(fill = "grey60", size = 4.5))) +
  scale_x_reverse(breaks = BR, expand = expansion(add = 2)) +
  scale_y_reverse(breaks = BR, expand = expansion(add = 2)) +
  labs(x = "Rank (no discordant-mismatch filter)", y = "Rank (Main analysis)") +
  theme_manuscript(base_size = 13, facets = TRUE) + grid_both +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(2, "pt"))
save_figure_sd(p, STEM, width = 9.5, height = 7.5, parts = c("ranks", "top_k", "rho"), grid = "both", script = SCRIPT)
