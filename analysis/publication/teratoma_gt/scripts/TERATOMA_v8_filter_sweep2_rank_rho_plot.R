#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_filter_sweep2_rank_rho_plot.R   (2026-10-21; draws from source data only since 2026-10-25)
# Response figure 6: sensitivity of the configuration ranking to the downstream filters, t = 50% KO truths.
#   top     Spearman rho of the overall configuration ranking vs the ranking under the main filters
#   bottom  rank, under the moved filter, of the configuration ranked first under the main filters
# Input : source_data/pm_gt/TERATOMA_v8_filter_sweep2_rank_rho_pwcov.t50.tsv
#         (TERATOMA_v8_filter_sweep2_rank_rho_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_filter_sweep2_rank_rho_pwcov.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_filter_sweep2_rank_rho_plot.R"
STEM <- "TERATOMA_v8_filter_sweep2_rank_rho_pwcov"

d <- read_part(STEM, "t50")
GATE_COL <- setNames(c("#0072B2", "#009E73"), levels(d$Gate))
relab <- function(x) sub("^[a-z]+ ", "", x)

row <- function(p, ylab, x_text)
  p + facet_grid(. ~ Axis, scales = "free_x", space = "free_x", labeller = label_value) +
    scale_x_discrete(labels = relab, expand = expansion(add = 0.75)) + scale_colour_manual(values = GATE_COL, name = NULL) +
    labs(x = NULL, y = ylab) + theme_manuscript(base_size = 14, facets = TRUE) + grid_y +
    theme(panel.spacing.x = unit(0.9, "lines")) +
    (if (!x_text) theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else theme())
p_top <- row(ggplot(d, aes(x, rho, colour = Gate, group = Gate)) +
               geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
               geom_point(data = d[centre == TRUE], shape = 21, size = 3.8, stroke = 1.1, fill = NA, colour = "black") +
               geom_text(data = d[centre == FALSE & !is.na(rho)], aes(label = sprintf("%.3f", rho), vjust = vj),
                         size = 4.3, show.legend = FALSE) +
               scale_y_continuous(limits = c(min(0, floor(min(d$rho, na.rm = TRUE) * 4) / 4), 1.08), breaks = seq(-1, 1, 0.25)),
             "Spearman's \u03C1 with the\nmain-filter ranking", x_text = FALSE) +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1))
YMAX <- max(3, max(d$centre_top_rank, na.rm = TRUE))
p_bot <- row(ggplot(d, aes(x, centre_top_rank, colour = Gate, group = Gate)) +
               geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
               geom_point(data = d[centre == TRUE], shape = 21, size = 3.8, stroke = 1.1, fill = NA, colour = "black") +
               scale_y_reverse(limits = c(YMAX + 0.3, 0.7), breaks = if (YMAX <= 6) 1:YMAX else unique(c(1, pretty(c(1, YMAX))[-1]))),
             "Rank of the main-filter\ntop configuration", x_text = TRUE) +
  theme(legend.position = "none", strip.text.x = element_blank(), strip.background = element_blank())
p <- (p_top / p_bot) + plot_layout(heights = c(1.2, 1))
save_figure_sd(p, STEM, width = 13.5, height = 7, parts = c("t50", "all_t"), grid = "y", script = SCRIPT)
