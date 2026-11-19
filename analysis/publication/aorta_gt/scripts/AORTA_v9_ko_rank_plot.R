#!/usr/bin/env Rscript
# ============================================================================
# AORTA_v9_ko_rank_plot.R   (2026-11-19)
#   a  aorta rank of the 76 configurations under each KO truth set (layout of response figure 7a)
#   b  teratoma vs aorta mean F1 across the six KO truth sets (least-squares line; shapes as response figure 13)
# Input : source_data/aorta_gt/AORTA_v9_ko_rank.{ranks,configs,cross}.tsv (AORTA_v9_ko_rank_source_data.R)
# Output: figures/aorta_gt/AORTA_v9_ko_rank.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "AORTA_plot_common.R"))
suppressPackageStartupMessages({ library(patchwork); library(ggrepel) })
SCRIPT <- "publication/aorta_gt/scripts/AORTA_v9_ko_rank_plot.R"
STEM <- "AORTA_v9_ko_rank"
BASE <- 12

rk <- read_part(STEM, "ranks")
cfg <- read_part(STEM, "configs")
cr <- read_part(STEM, "cross")
N <- nrow(cfg)
TRUTH_LAB <- V8_TRUTH_LABELS[V8_KO_TRUTHS]
cmap_truth <- setNames(unname(V8_TRUTH_COLOURS[V8_KO_TRUTHS]), unname(TRUTH_LAB))
rk[, Truth := factor(Truth, TRUTH_LAB)]
CMAPS <- list(Preprocessing = cmap_proc, Caller = cmap_method, Stranding = cmap_strand, `SNP removal` = cmap_snp)
for (tb in list(cfg, cr)) {
  tb[, `SNP removal` := relabel(snp_subtract, label_snpSubtract)]
  stopifnot(all(tb$Preprocessing %in% names(cmap_proc)), all(tb$Caller %in% names(cmap_method)),
            all(tb$Stranding %in% names(cmap_strand)), all(tb$`SNP removal` %in% names(cmap_snp)))
}
CALLERS <- sort(unique(cfg$Caller))

# ---- a: aorta ranks ---------------------------------------------------------------------
Y_LIM <- c(0.5, N + 0.5)
tile_col <- function(col) {
  d <- cfg[, .(y, value = get(col), column = col)]
  brk <- if (col == "Caller") CALLERS else intersect(names(CMAPS[[col]]), unique(d$value))
  ggplot(d, aes(column, y, fill = value)) +
    geom_tile(width = 1, height = 0.78) +
    scale_fill_manual(values = CMAPS[[col]], breaks = brk, name = col) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 1, title.position = "top", order = match(col, names(CMAPS)))) +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = BASE) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
          axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(5, 0, 5, 0))
}
d <- copy(rk)[, y := y + (match(truth, V8_KO_TRUTHS) - 3.5) * 0.06]
p_dots <- ggplot(d, aes(overall_rank, y)) +
  geom_segment(data = cfg, aes(x = best, xend = worst, y = y, yend = y), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = Truth), size = 2) +
  scale_x_continuous(breaks = c(1, seq(5, N, 5)), limits = c(1, N + 1), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
  scale_colour_manual_cvd(values = cmap_truth, name = "Truth set") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", order = 5, override.aes = list(size = 3))) +
  labs(x = "Overall cross-tier rank (aorta)", y = NULL) +
  theme_manuscript(base_size = BASE) + grid_x +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey80"), plot.margin = margin(5, 5, 5, 4))
TILE_IN <- 0.22
p_a <- wrap_plots(c(lapply(names(CMAPS), tile_col), list(p_dots)), nrow = 1,
                  widths = unit(c(rep(TILE_IN, length(CMAPS)), 1), c(rep("in", length(CMAPS)), "null")))

# ---- c: teratoma vs aorta (point style of response figure 13) ------------------------------
shp_proc <- setNames(c(21, 22, 23, 24, 25, 4),
                     c("De-duplicated (coordinate)", "De-duplicated (umi_tools)", "De-duplicated (UMI-aware MarkDuplicates)",
                       "Recalibrated", "Recalibrated (Q20)", "Recalibrated (default)"))
stopifnot(setequal(names(shp_proc), cr$Preprocessing))
cr[, Preprocessing := factor(Preprocessing, names(shp_proc))]
p_c <- ggplot(cr, aes(mean_f1_teratoma, mean_f1_aorta)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey30", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(fill = Caller, shape = Preprocessing), size = 2.6, stroke = 0.3, colour = "grey15", alpha = 0.9) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5, size = 4.2,
           label = paste0(spearman_lab(cr$rho[1]), "\nn = ", N, " configurations")) +
  scale_fill_manual(values = cmap_method, breaks = CALLERS, guide = "none") +
  scale_shape_manual(values = shp_proc, guide = "none") +
  labs(x = "Mean F1 across six KO truth sets\n(human teratoma)", y = "Mean F1 across six KO truth sets\n(mouse aorta)") +
  theme_manuscript(base_size = BASE) +
  theme(panel.grid.major = element_line(linewidth = 0.3, colour = "grey88"), aspect.ratio = 1)

p <- wrap_plots(p_a, wrap_plots(p_c, plot_spacer(), ncol = 1, heights = c(1, 0.45)), nrow = 1, widths = c(1, 1.1)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("a", rep("", length(CMAPS)), "b"))) &
  theme(legend.position = "bottom", legend.justification = "left", legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = BASE), legend.text = element_text(size = BASE - 1),
        plot.tag = element_text(face = "bold", size = BASE + 6))
save_figure_sd(p, STEM, width = 14, height = 12, parts = c("ranks", "configs", "cross"),
               grid = "x", script = SCRIPT)
