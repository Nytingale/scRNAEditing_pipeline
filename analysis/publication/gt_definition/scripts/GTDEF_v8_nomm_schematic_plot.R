#!/usr/bin/env Rscript
# =============================================================================
# GTDEF_v8_nomm_schematic_plot.R   (2026-10-23; Reviewer 2, major 3: what the removal of sites near other mismatch
# types does; redrawn from original manuscript Fig. 1e, window 1)
# A site is removed only when a mismatch of a different type lies within +- w nt on the same strand. Neighbouring
# sites of the same type (a cluster of A-to-G sites, as in hyper-editing) never trigger removal; isolated sites are
# also kept. Positions are illustrative and set below, so this figure has no source data table.
# Output: publication/gt_definition/figures/GTDEF_v8_nomm_schematic.{svg,png},
#         review/response_figures_2026-10-15/GTDEF_v8_nomm_schematic.png (response figure 14)
# =============================================================================
B <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "..", "..", ".."))
source(file.path(B, "plotting", "plotting_functions.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
FIG <- file.path(B, "publication/gt_definition/figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
W <- 100
X_MAX <- 1020

s <- data.table(pos = c(60, 95, 125, 160, 190, 520, 575, 900),
                genome = c("A", "A", "A", "A", "A", "A", "C", "A"),
                rna = c("G", "G", "G", "G", "G", "G", "T", "G"),
                group = c(rep("cluster", 5), "pair", "pair", "isolated"))
s[, type := paste0(genome, "-to-", rna)]
s[, near_other := vapply(seq_len(.N), function(i) any(abs(pos[-i] - pos[i]) <= W & type[-i] != type[i]), logical(1))]
s[, outcome := factor(fifelse(near_other, "Removed", "Kept"), names(cmap_filter_keep))]

win <- data.table(pos = c(125, 520, 575, 900), y = c(3.0, 3.0, 3.45, 3.0))
win[, `:=`(lo = pmax(0, pos - W), hi = pmin(X_MAX, pos + W))]
grp <- data.table(x = c(125, 547.5, 900), y = -0.75,
                  lab = c("Cluster of A-to-G sites\n(e.g. hyper-editing): all kept",
                          "A-to-G site with a different\nmismatch type within \u00b1 w nt: both removed",
                          "Isolated A-to-G site:\nkept"))
rows <- data.table(y = c(2, 1), lab = c("Genome", "RNA"))

p <- ggplot() +
  geom_segment(data = rows, aes(x = 0, xend = X_MAX, y = y, yend = y), linewidth = 1.1) +
  geom_segment(data = s, aes(x = pos, xend = pos, y = 1.25, yend = 3.4), linetype = "dashed", colour = "grey60",
               linewidth = 0.4) +
  geom_segment(data = win, aes(x = lo, xend = hi, y = y, yend = y), colour = okabe_ito[["skyblue"]], linewidth = 1.2,
               arrow = arrow(ends = "both", length = unit(0.12, "inches"), type = "closed")) +
  annotate("text", x = 125, y = 3.32, label = "\u00b1 w nt", size = 5.5) +
  geom_label(data = s, aes(pos, 2, label = genome), size = 7, linewidth = 0, fill = "white",
             label.padding = unit(0.08, "lines")) +
  geom_label(data = s, aes(pos, 1, label = rna), size = 7, linewidth = 0, fill = "white",
             label.padding = unit(0.08, "lines")) +
  geom_point(data = s, aes(pos, 0.15, colour = outcome, shape = outcome), size = 6, stroke = 2.2) +
  geom_text(data = grp, aes(x, y, label = lab), size = 5.2, lineheight = 0.95) +
  geom_text(data = rows, aes(x = -15, y = y, label = lab), hjust = 1, size = 6) +
  scale_colour_manual_cvd(values = cmap_filter_keep, name = NULL) +
  scale_shape_manual(values = c(Kept = 16, Removed = 4), name = NULL) +
  scale_x_continuous(limits = c(-120, X_MAX + 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.3, 3.7), expand = c(0, 0)) +
  theme_manuscript(base_size = 16) + theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 16), plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(FIG, "GTDEF_v8_nomm_schematic.svg"), p, width = 13, height = 4.6, device = grDevices::svg)
ggsave(file.path(FIG, "GTDEF_v8_nomm_schematic.png"), p, width = 13, height = 4.6, dpi = 200)
dir.create(file.path(B, "review/response_figures_2026-10-15"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(B, "review/response_figures_2026-10-15/GTDEF_v8_nomm_schematic.png"), p, width = 13, height = 4.6,
       dpi = 200)
print(s)
