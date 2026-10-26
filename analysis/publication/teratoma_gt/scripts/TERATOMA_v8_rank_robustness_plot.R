#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_rank_robustness_plot.R   (2026-10-13; draws from source data only since 2026-10-25)
# Rank of every configuration under the nine v8 truth sets (six KO ground truths, three REDIportal).
#   all_tiles  response figure 4: configuration tiles and the rank under each truth set (+ _legend)
#   by_family  response figure 5: a within the KO ground truths, b within REDIportal
# Input : source_data/pm_gt/TERATOMA_v8_rank_robustness.{ranks,configs,by_family}.tsv
#         (TERATOMA_v8_rank_robustness_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_rank_robustness_all_tiles{,_legend}.svg
#         figures/pm_gt/TERATOMA_v8_rank_robustness_by_family.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_rank_robustness_plot.R"
STEM <- "TERATOMA_v8_rank_robustness"

rk <- read_part(STEM, "ranks")
all_cfg <- read_part(STEM, "configs")
fam_rank <- read_part(STEM, "by_family")

COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
cmap_truth <- setNames(unname(V8_TRUTH_COLOURS[V8_TRUTH_IDS]), unname(V8_TRUTH_LABELS[V8_TRUTH_IDS]))
N_TRUTH <- length(V8_TRUTH_IDS)
N <- nrow(all_cfg)

# ---- tiles figure ----------------------------------------------------------------
CMAPS <- list(Preprocessing = cmap_proc, Caller = cmap_method, Stranding = cmap_strand, `SNP removal` = cmap_snp)
Y_LIM_ALL <- c(0.5, N + 0.5)
BASE_ALL <- 20
tile_col <- function(col) {
  d <- all_cfg[, .(y, value = get(col), column = col)]
  ggplot(d, aes(column, y, fill = value)) +
    geom_tile(width = 1, height = 0.78) +
    scale_fill_manual(values = CMAPS[[col]], name = col,
                      breaks = if (col == "Caller") {
                        v <- unique(d$value); v[order(tolower(v))]
                      } else intersect(names(CMAPS[[col]]), unique(d$value))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 1, title.position = "top")) +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = BASE_ALL) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
          axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(5, 0, 5, 0))
}
d <- copy(rk)
d[, y := (N + 1L - match(as.character(label), as.character(all_cfg$label))) +
        (as.integer(Truth) - (N_TRUTH + 1) / 2) * 0.04]
cf <- copy(all_cfg)
p_dots_all <- ggplot(d, aes(overall_rank, y)) +
  geom_segment(data = cf, aes(x = best, xend = worst, y = y, yend = y), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = Truth), size = 4.6) +
  geom_text(data = cf, aes(x = worst, y = y, label = worst), nudge_x = 0.9, hjust = 0, size = 7, fontface = "bold") +
  scale_x_continuous(breaks = c(1, seq(5, 55, 5)), expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
  scale_colour_manual_cvd(values = cmap_truth, name = NULL) +
  guides(colour = guide_legend(ncol = 1, title = "Truth set", title.position = "top", override.aes = list(size = 5))) +
  labs(x = "Overall cross-tier rank", y = NULL) +
  theme_manuscript(base_size = BASE_ALL) + grid_x +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey77"), plot.margin = margin(5, 5, 5, 4))
TILE_IN <- 0.36
legend_theme <- function(text_size, key_cm, just = "top")
  theme(legend.position = "bottom", legend.justification = just,
        legend.title = element_text(face = "bold", size = text_size + 1),
        legend.text = element_text(size = text_size), legend.key.size = unit(key_cm, "cm"))
panels <- c(lapply(names(CMAPS), tile_col), list(p_dots_all))
p_tiles <- patchwork::wrap_plots(panels, nrow = 1,
  widths = unit(c(rep(TILE_IN, length(CMAPS)), 1), c(rep("in", length(CMAPS)), "null"))) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = paste0("Rank of ", N, " method configurations under ", N_TRUTH, " truth definitions"),
    subtitle = V8_SUBTITLE,
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
                  plot.subtitle = element_text(hjust = 0.5, size = 14))) &
  legend_theme(10, 0.35)

legend_stack <- function(text_size, key_cm, pt_size) {
  legs <- lapply(seq_along(panels), function(i) {
    pp <- panels[[i]] + legend_theme(text_size, key_cm, just = "left")
    if (i == length(panels))
      pp <- pp + guides(colour = guide_legend(ncol = 1, title = "Truth set", title.position = "top",
                                              override.aes = list(size = pt_size)))
    cowplot::get_plot_component(pp, "guide-box-bottom")
  })
  cowplot::plot_grid(plotlist = legs, ncol = 1, rel_heights = leg_h, align = "v", axis = "l")
}
n_keys <- c(sapply(names(CMAPS), function(col) uniqueN(all_cfg[[col]])), N_TRUTH)
leg_h <- n_keys + 2.2
p_leg <- legend_stack(18, 0.8, 7)
save_figure_sd(p_leg, "TERATOMA_v8_rank_robustness_all_tiles_legend", width = 11, height = 0.4 * sum(leg_h),
               parts = "configs", grid = "none", script = SCRIPT, src_stem = STEM)
ALL_H <- 3 + TILE_IN * N
LEG_H <- 0.52 * sum(leg_h)
PAD_H <- max(ALL_H - LEG_H, 0) / 2
p_all <- cowplot::plot_grid(
  patchwork::patchworkGrob(p_tiles & theme(legend.position = "none", plot.title = element_blank(),
                                            plot.subtitle = element_blank())),
  cowplot::plot_grid(NULL, legend_stack(24, 1.05, 9), NULL, ncol = 1, rel_heights = c(PAD_H, LEG_H, PAD_H)),
  nrow = 1, rel_widths = c(16, 14))
save_figure_sd(p_all, "TERATOMA_v8_rank_robustness_all_tiles", width = 30, height = ALL_H,
               parts = c("ranks", "configs"), grid = "x", script = SCRIPT, src_stem = STEM)

# ---- the same plot within each truth family ----------------------------------------
FAMILIES <- list(ko = list(ids = V8_KO_TRUTHS, title = "Editing reduction in KO ground truths"),
                 redi = list(ids = paste0("redi_", V8_REDI_VARIANTS), title = "REDIportal ground truths"))
fam_tile <- function(cfg, col, legend = FALSE) {
  d <- cfg[, .(y, value = get(col), column = col)]
  ggplot(d, aes(column, y, fill = value)) +
    geom_tile(width = 1, height = 0.78) +
    scale_fill_manual(values = CMAPS[[col]], name = col, guide = if (legend) "legend" else "none",
                      breaks = if (col == "Caller") { v <- unique(d$value); v[order(tolower(v))] }
                               else intersect(names(CMAPS[[col]]), unique(d$value))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
    guides(fill = if (legend) guide_legend(ncol = 1, title.position = "top") else "none") +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = BASE_ALL) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
          axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(5, 0, 5, 0))
}
fam_panel <- function(nm) {
  fm <- FAMILIES[[nm]]
  ids <- fm$ids
  cfg <- fam_rank[fam == nm]
  dd <- merge(rk[truth %in% ids, c(COMBO, "truth", "Truth", "overall_rank"), with = FALSE],
              cfg[, c(COMBO, "y"), with = FALSE], by = COMBO)
  dd[, `:=`(Truth = droplevels(Truth), y = y + (match(truth, ids) - (length(ids) + 1) / 2) * 0.06)]
  p_dots <- ggplot(dd, aes(overall_rank, y)) +
    geom_segment(data = cfg, aes(x = best, xend = worst, y = y, yend = y), colour = "grey70", linewidth = 0.5) +
    geom_point(aes(colour = Truth), size = 4.6) +
    geom_text(data = cfg, aes(x = worst, y = y, label = worst), nudge_x = 0.9, hjust = 0, size = 7, fontface = "bold") +
    scale_x_continuous(breaks = c(1, seq(5, 55, 5)), limits = c(1, N + 2), expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
    scale_colour_manual_cvd(values = cmap_truth[levels(dd$Truth)], name = "Truth set") +
    guides(colour = "none") +
    labs(x = "Overall cross-tier rank", y = NULL) +
    theme_manuscript(base_size = BASE_ALL) + grid_x +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey77"), plot.margin = margin(5, 5, 5, 4))
  patchwork::wrap_plots(c(lapply(names(CMAPS), function(col) fam_tile(cfg, col)), list(p_dots)), nrow = 1,
         widths = unit(c(rep(TILE_IN, length(CMAPS)), 1), c(rep("in", length(CMAPS)), "null"))) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = fm$title,
                               theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
                                             plot.margin = margin(40, 10, 10, 10)))
}
fp <- setNames(lapply(names(FAMILIES), fam_panel), names(FAMILIES))
p_fam <- cowplot::plot_grid(patchwork::patchworkGrob(fp$ko), patchwork::patchworkGrob(fp$redi),
                            nrow = 1, labels = c("a", "b"), label_size = 44)
FAM_H <- 3 + TILE_IN * N
FAM_PAD <- max(FAM_H - LEG_H, 0) / 2
p_fam <- cowplot::plot_grid(p_fam,
  cowplot::plot_grid(NULL, legend_stack(24, 1.05, 9), NULL, ncol = 1, rel_heights = c(FAM_PAD, LEG_H, FAM_PAD)),
  nrow = 1, rel_widths = c(28, 14))
save_figure_sd(p_fam, "TERATOMA_v8_rank_robustness_by_family", width = 42, height = FAM_H,
               parts = c("ranks", "by_family"), grid = "x", script = SCRIPT, src_stem = STEM)
