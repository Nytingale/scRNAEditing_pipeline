#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_rank_robustness_plot.R   (HANDOFF 2026-09-29 UPDATE-26; nine truths 2026-10-01)
# Overall cross-tier rank of each reported configuration (full benchmark call set) under the nine
# truth definitions in TERATOMA_common.R (NINE_TRUTHS), all global (one list of true sites for
# every unit); configurations ordered by their worst rank across the nine.
#   ground truth at 25 / 50 / 75% depletion, KO arm covered if expected KO alt >= 3
#   ground truth at 25 / 50 / 75% depletion, KO arm covered if NTC and KO cov >= 10
#   REDIportal (not dbSNP common), merged NTC cov >= 10 with >= 1 G / >= 3 G / >= 3 G and edit >= 5%
# Rank = F1 rank within tier summed over tiers, ties broken by higher mean F1 (as in the grids).
#
# Inputs
#   human/pm/gt_eval_v6/scoring/sweep_infneg{,_depth}/teratoma_inputs/invariance/d0.{25,50,75}/gt_tier_scoring.tsv
#   human/pm/gt_eval_v6/rediportal/scoring/gt_tier_scoring_{cov10_alt1,alt3,alt3_e05}.tsv
# Outputs (figures/pm_gt/)
#   TERATOMA_rank_robustness.svg        top configurations by worst rank
#   TERATOMA_rank_robustness_all.svg    all reported configurations (supplementary)
#   TERATOMA_rank_robustness_all_tiles.svg (+ _legend.svg)  all configurations, colour tiles; legends separate
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_rank_robustness_plot.R"

COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
TRUTH_LAB <- setNames(NINE_TRUTHS$label, NINE_TRUTHS$truth)
cmap_truth <- setNames(NINE_TRUTHS$colour, NINE_TRUTHS$label)
SRC <- setNames(NINE_TRUTHS$file, NINE_TRUTHS$truth)
N_TRUTH <- length(TRUTH_LAB)

overall_rank <- function(s) {
  s <- copy(s)[, f1 := fifelse(precision + recall_vs_union > 0,
                               2 * precision * recall_vs_union / (precision + recall_vs_union), NA_real_)]
  s[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
  r <- s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
  r[, overall_rank := frank(list(rank_sum, -mean_f1), ties.method = "min")]
}
rk <- rbindlist(lapply(names(SRC), function(n) overall_rank(read_truth_scoring(n))[, truth := n]))
rk[, Truth := factor(TRUTH_LAB[truth], TRUTH_LAB)]

cfg <- rk[, .(worst = max(overall_rank), best = min(overall_rank), mean_rank = mean(overall_rank)), by = COMBO]
setorder(cfg, worst, mean_rank)
cfg[, label := sprintf("%s | %s | %s | %s", relabel(caller, label_caller), relabel(bam_type, label_bamType),
                       relabel(strand_mode, label_strandMode), relabel(snp_subtract, label_snpSubtract))]
cfg[, label := factor(label, rev(label))]
rk <- merge(rk, cfg, by = COMBO)
fwrite(rk[order(worst, mean_rank, Truth)], file.path(SRC_DIR, "TERATOMA_rank_robustness.tsv"), sep = "\t")

rank_plot <- function(k, jitter = 0.08, base_size = 11, pt_size = 2.4, lab_size = 3.6, x_breaks = c(1, 5, 10, 15, 20, 30, 40, 50)) {
  keep <- cfg$label[seq_len(min(k, nrow(cfg)))]
  d <- rk[label %in% keep]; cf <- cfg[label %in% keep]
  d[, y := as.integer(droplevels(label)) + (as.integer(Truth) - (N_TRUTH + 1) / 2) * jitter]
  cf[, y := as.integer(droplevels(label))]
  ggplot(d, aes(overall_rank, y)) +
    geom_segment(data = cf, aes(x = best, xend = worst, y = y, yend = y), colour = "grey70", linewidth = 0.5) +
    geom_point(aes(colour = Truth), size = pt_size) +
    geom_text(data = cf, aes(x = worst, y = y, label = worst), nudge_x = 0.9, hjust = 0, size = lab_size, fontface = "bold") +
    scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0.02, 0.06))) +
    scale_y_continuous(breaks = seq_along(levels(droplevels(cf$label))), labels = levels(droplevels(cf$label)),
                       expand = expansion(add = 0.5)) +
    scale_colour_manual_cvd(values = cmap_truth, name = NULL) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE)) +
    labs(x = "Overall cross-tier rank", y = NULL,
         title = paste0("Rank of each configuration under ", N_TRUTH, " truth definitions")) +
    theme_manuscript(base_size = base_size) + grid_x + theme(legend.location = "plot", plot.title.position = "plot")
}

# Top-K figure: configuration shown as the cross-tier grid's coloured tiles (same palettes and
# text-colour rule), dot panel on the same row positions.
TOP_K <- 15L
top <- cfg[seq_len(TOP_K)][, y := TOP_K + 1L - .I]
GRID_COLS <- c("Preprocessing", "Caller", "Stranding", "SNP removal")
top[, `:=`(Preprocessing = relabel(bam_type, label_bamType), Caller = relabel(caller, label_caller),
           Stranding = relabel(strand_mode, label_strandMode), `SNP removal` = relabel(snp_subtract, label_snpSubtract))]
tiles <- melt(top[, c("y", GRID_COLS), with = FALSE], id.vars = "y", variable.name = "column", value.name = "value")
tiles[, column := factor(column, GRID_COLS)]
tiles[, fill := fcase(column == "Preprocessing", cmap_proc[value], column == "Caller", cmap_method[value],
                      column == "Stranding", cmap_strand[value], column == "SNP removal", cmap_snp[value])]
tiles[, txt_col := {rgb <- col2rgb(fill); ifelse(rgb[1, ] * 0.299 + rgb[2, ] * 0.587 + rgb[3, ] * 0.114 < 120, "white", "black")}]
Y_LIM <- c(0.5, TOP_K + 0.5)
p_tiles <- ggplot(tiles, aes(column, y)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.3) +
  geom_text(aes(label = value, colour = txt_col), size = 3.6) +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_discrete(position = "top") +
  scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, face = "bold"), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = margin(5, 0, 5, 5))
p_dots <- rank_plot(TOP_K) +
  scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
  labs(title = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "none", plot.margin = margin(5, 5, 5, 0))
legend <- cowplot::get_plot_component(rank_plot(TOP_K), "guide-box-bottom", return_all = TRUE)
if (is.list(legend) && !inherits(legend, "grob")) legend <- legend[[1]]
title_grob <- cowplot::ggdraw() +
  cowplot::draw_label(paste0("Rank of each configuration under ", N_TRUTH, " truth definitions"), fontface = "bold", size = 13)
body <- cowplot::plot_grid(p_tiles, p_dots, ncol = 2, rel_widths = c(1.5, 1), align = "h", axis = "tb")
p <- cowplot::plot_grid(title_grob, body, legend, ncol = 1, rel_heights = c(0.4, 0.42 * TOP_K + 0.9, 0.8))
save_figure(p, "TERATOMA_rank_robustness", width = 16, height = 1.2 + 0.42 * TOP_K + 0.9,
            sources = SRC, data = rk[label %in% cfg$label[1:TOP_K]], grid = "x", script = SCRIPT)
p <- rank_plot(nrow(cfg))
save_figure(p, "TERATOMA_rank_robustness_all", width = 13, height = 1.8 + 0.36 * nrow(cfg),
            sources = SRC, data = rk, grid = "x", script = SCRIPT)

# All configurations, configuration as unlabelled square tiles (one legend per column) ----
N <- nrow(cfg)
all_cfg <- copy(cfg)[, y := N + 1L - .I]
all_cfg[, `:=`(Preprocessing = relabel(bam_type, label_bamType), Caller = relabel(caller, label_caller),
               Stranding = relabel(strand_mode, label_strandMode), `SNP removal` = relabel(snp_subtract, label_snpSubtract))]
CMAPS <- list(Preprocessing = cmap_proc, Caller = cmap_method, Stranding = cmap_strand, `SNP removal` = cmap_snp)
Y_LIM_ALL <- c(0.5, N + 0.5)
BASE_ALL <- 20
tile_col <- function(col) {
  d <- all_cfg[, .(y, value = get(col), column = col)]
  ggplot(d, aes(column, y, fill = value)) +
    geom_tile(width = 1, height = 0.78) +
    scale_fill_manual(values = CMAPS[[col]], name = col,
                      breaks = if (col == "Caller") {v <- unique(d$value); v[order(tolower(v))]}
                               else intersect(names(CMAPS[[col]]), unique(d$value))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 1, title.position = "top")) +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = BASE_ALL) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
          axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(5, 0, 5, 0))
}
p_dots_all <- rank_plot(N, jitter = 0.04, base_size = BASE_ALL, pt_size = 4.6, lab_size = 7, x_breaks = c(1, seq(5, 55, 5))) +
  scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
  labs(title = NULL) +
  guides(colour = guide_legend(ncol = 1, title = "Truth set", title.position = "top",
                               override.aes = list(size = 5))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey82"),
        plot.margin = margin(5, 5, 5, 4))
TILE_IN <- 0.36
legend_theme <- function(text_size, key_cm, just = "top")
  theme(legend.position = "bottom", legend.justification = just,
        legend.title = element_text(face = "bold", size = text_size + 1), legend.text = element_text(size = text_size),
        legend.key.size = unit(key_cm, "cm"))
panels <- c(lapply(names(CMAPS), tile_col), list(p_dots_all))

# main figure: small legends, large axes and points
p <- patchwork::wrap_plots(panels, nrow = 1,
                           widths = unit(c(rep(TILE_IN, length(CMAPS)), 1), c(rep("in", length(CMAPS)), "null"))) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = paste0("Rank of ", N, " method configurations under ", length(TRUTH_LAB), " truth definitions"),
                             theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22))) &
  legend_theme(10, 0.35)
save_figure(p, "TERATOMA_rank_robustness_all_tiles", width = 16, height = 2.4 + TILE_IN * N + 2.2,
            sources = SRC, data = all_cfg, grid = "x", script = SCRIPT)

# legend-only figure (large keys and text)
legs <- lapply(seq_along(panels), function(i) {
  pp <- panels[[i]] + legend_theme(18, 0.8, just = "left")
  if (i == length(panels)) pp <- pp + guides(colour = guide_legend(ncol = 1, title = "Truth set", title.position = "top",
                                                                   override.aes = list(size = 7)))
  cowplot::get_plot_component(pp, "guide-box-bottom")
})
n_keys <- c(sapply(names(CMAPS), function(col) uniqueN(all_cfg[[col]])), length(TRUTH_LAB))
leg_h <- n_keys + 2.2
p_leg <- cowplot::plot_grid(plotlist = legs, ncol = 1, rel_heights = leg_h, align = "v", axis = "l")
save_figure(p_leg, "TERATOMA_rank_robustness_all_tiles_legend", width = 8, height = 0.4 * sum(leg_h),
            sources = SRC, data = all_cfg, grid = "none", script = SCRIPT)

print(cfg[1:TOP_K, .(label, best, worst, mean_rank = round(mean_rank, 1))])
