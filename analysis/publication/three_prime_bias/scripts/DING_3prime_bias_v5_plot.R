#!/usr/bin/env Rscript
# =============================================================================
# DING_3prime_bias_v5_plot.R   (2026-10-24; split from human/scripts/ding_3prime_bias_fig_v5.R 2026-10-25;
# source data only)
# Response figure 16 (Reviewer 3, major 2): Ding et al. human PBMC1, 10x Chromium 3' v3 vs Smart-seq2 at the same
# number of unique reads. One pair of panels per gene region: left, edited sites per platform (Smart-seq2 = mean over
# 3 subsampling seeds, bar = min-max); right, edited sites by distance from the gene 3' end as % of each platform's
# edited sites in the region. b: edited sites found by Smart-seq2 only, both platforms or 10x only, per region of the
# site annotation (area-proportional; mean over the 3 seeds; user 2026-10-25, as panel b of v2).
# Input : publication/three_prime_bias/source_data/ding_3prime_bias_v5.{edited_by_region,dist3p_bins,overlap_by_region}.tsv
#         (DING_3prime_bias_v5_source_data.R)
# Output: publication/three_prime_bias/figures/ding_3prime_bias_v5.{svg,png},
#         review/response_figures_2026-10-15/DING_3prime_bias_v5.png
# =============================================================================
B <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "..", "..", ".."))
source(file.path(B, "plotting", "plotting_functions.R"))
source(file.path(B, "plotting", "source_data_io.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
PUB <- file.path(B, "publication/three_prime_bias")
SRC_DIR <- file.path(PUB, "source_data")
STEM <- "ding_3prime_bias_v5"
for (d in c(file.path(PUB, "figures"), file.path(B, "review/response_figures_2026-10-15")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

cnt <- read_sd(SRC_DIR, STEM, "edited_by_region")
c_dt <- read_sd(SRC_DIR, STEM, "dist3p_bins")
GB_REG <- levels(cnt$region)

pair <- function(r, first_col, bottom_row) {
  pa <- ggplot(cnt[region == r], aes(platform, N, fill = platform)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = N_min, ymax = N_max), width = 0.25, linewidth = 0.5) +
    geom_text(aes(label = scales::comma(round(N)), y = N_max), vjust = -0.9, size = 5) +
    scale_y_continuous(labels = shorten_axis, expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = cmap_3prime_platform, guide = "none") +
    labs(x = NULL, y = if (first_col) "Edited sites" else NULL) +
    theme_manuscript(base_size = 16) + grid_y +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4, 4, 4, 12))
  pc <- ggplot(c_dt[region == r], aes(x, pct, colour = platform, fill = platform)) +
    geom_ribbon(aes(ymin = pct_min, ymax = pct_max), colour = NA, alpha = 0.3) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    scale_x_log10(breaks = 10^c(0, 2, 4, 6), labels = shorten_axis, expand = expansion(mult = 0.03)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.08))) +
    scale_colour_manual_cvd(values = cmap_3prime_platform, name = NULL) +
    scale_fill_manual(values = cmap_3prime_platform, guide = "none") +
    labs(x = if (bottom_row) "Distance from gene 3' end (nt)" else NULL,
         y = if (first_col) "Edited sites in the region (%)" else NULL) +
    theme_manuscript(base_size = 16) + grid_both +
    theme(plot.margin = margin(4, 30, 4, 4))
  # one strip across both panels, so the pair reads as one facet
  strip <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = r, size = 6.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(14, 30, 2, 12))
  strip / ((pa | pc) + plot_layout(widths = c(1, 2.6))) + plot_layout(heights = c(0.09, 1))
}
NCOL <- 2
pairs <- lapply(seq_along(GB_REG), function(i) {
  pair(GB_REG[i], first_col = (i - 1) %% NCOL == 0, bottom_row = i > length(GB_REG) - NCOL)
})
cells <- unlist(lapply(seq(1, length(pairs), by = NCOL), function(i) list(pairs[[i]], plot_spacer(), pairs[[i + 1]])),
                recursive = FALSE)
fig_a <- wrap_plots(cells, ncol = 2 * NCOL - 1, widths = c(1, 0.04, 1)) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.key.width = unit(2.5, "lines"))

# b: area-proportional two-circle Euler diagram per region (each region on its own scale, union area 1); the circle
# distance is solved so that the lens area equals the shared sites
lens <- function(d, r1, r2) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1)) + r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2)) -
    0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
}
ov <- dcast(read_sd(SRC_DIR, STEM, "overlap_by_region"), region ~ set, value.var = "N")
geo <- ov[, {
  u <- ss2_only + shared + x10_only
  r1 <- sqrt((ss2_only + shared) / u / pi); r2 <- sqrt((x10_only + shared) / u / pi)
  d <- if (shared > 0) uniroot(function(z) lens(z, r1, r2) - shared / u, c(abs(r1 - r2) + 1e-9, r1 + r2 - 1e-9))$root
       else r1 + r2 + 0.05
  .(r1 = r1, r2 = r2, d = d, ss2_only = ss2_only, shared = shared, x10_only = x10_only)
}, by = region]
PLAT <- names(cmap_3prime_platform)
th <- seq(0, 2 * pi, length.out = 200)
circ <- geo[, .(platform = rep(PLAT, each = length(th)),
                x = c(r1 * cos(th), d + r2 * cos(th)), y = c(r1 * sin(th), r2 * sin(th))), by = region]
circ[, platform := factor(platform, PLAT)]
lab <- geo[, {
  x_in <- c(mean(c(-r1, d - r2)), mean(c(d - r2, r1)), mean(c(r1, d + r2)))
  narrow <- (r1 - (d - r2)) < 0.3
  .(x = x_in,
    y = c(0, if (narrow) -max(r1, r2) - 0.12 else 0, 0),
    label = scales::comma(round(c(ss2_only, shared, x10_only))), narrow = narrow, xs = x_in[2],
    ys = 0)
}, by = region]
pb <- ggplot() +
  geom_polygon(data = circ, aes(x, y, fill = platform, colour = platform, group = platform), alpha = 0.45, show.legend = FALSE,
               linewidth = 0.5) +
  geom_segment(data = lab[narrow == TRUE & y < 0], aes(x = xs, xend = x, y = ys, yend = y + 0.05),
               linewidth = 0.3, colour = "grey30") +
  geom_text(data = lab, aes(x, y, label = label), size = 5) +
  # invisible points carry the legend, so its keys are circles rather than polygon squares
  geom_point(data = circ[, .SD[1], by = platform], aes(x, y, fill = platform), shape = 21, size = 0, alpha = 0) +
  facet_wrap(~ region, nrow = 2) +
  coord_fixed(clip = "off") +
  scale_fill_manual(values = cmap_3prime_platform, name = NULL,
                    guide = guide_legend(override.aes = list(size = 9, alpha = 0.45, colour = cmap_3prime_platform,
                                                             stroke = 0.6))) +
  scale_colour_manual(values = cmap_3prime_platform, guide = "none") +
  theme_void(base_size = 16) +
  theme(strip.text = element_text(face = "bold", size = 17, margin = margin(b = 6)),
        panel.spacing.y = unit(1.2, "lines"), plot.margin = margin(10, 10, 10, 10),
        legend.position = "bottom", legend.text = element_text(size = 18), legend.key.spacing.x = unit(1.5, "lines"))

fig <- wrap_elements(full = fig_a) / wrap_elements(full = pb) + plot_layout(heights = c(15, 7.3)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 22))
W <- 17; H <- 22.3
ggsave(file.path(PUB, "figures", "ding_3prime_bias_v5.svg"), fig, width = W, height = H, device = grDevices::svg)
ggsave(file.path(PUB, "figures", "ding_3prime_bias_v5.png"), fig, width = W, height = H, dpi = 200)
ggsave(file.path(B, "review/response_figures_2026-10-15", "DING_3prime_bias_v5.png"), fig, width = W, height = H, dpi = 200)
