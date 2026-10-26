#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_bam_processing_plot.R   (draws from source data only since 2026-10-25)
#
# Ground-truth DEFINITION panels on BAM processing (v6): mismatch types among sites called in NTC on each BAM type,
# before vs after restricting to the 12-type ground truth, and the overlap of ground-truth sites across BAM types.
# The coordinate vs UMI de-duplication hex panels are in GTDEF_hex_dedup_coord_vs_umi_{source_data,plot}.R.
#
# Input : source_data/pm_gt/GTDEF_bam_processing.<part>.tsv (GTDEF_bam_processing_source_data.R)
# Output: figures/pm_gt/GTDEF_bar_mismatch_spectrum{,6}_by_bamtype_vs_TRUE.svg,
#         GTDEF_euler_TRUE_sites_by_bamtype_<tier>.svg, GTDEF_upset_TRUE_sites_by_bamtype{,_<tier>}.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
suppressPackageStartupMessages({
  library(eulerr)
  library(ggh4x)
  library(egg)
  library(forcats)
})
SCRIPT <- "publication/gt_definition/scripts/GTDEF_bam_processing_plot.R"
STEM <- "GTDEF_bam_processing"
TRUTH_D <- 0.5
TIERS <- c("dataset", "sample", "celltype")
BT_REPORTED <- c("sampleRmDup", "sampleUmiDup", "recal", "recalq20", "recaldefault")
varType6_levels <- c("AC", "AG", "AT", "CA", "CG", "CT")

tier_title <- function(x) {
  x <- as.character(x)
  # Pipeline token `dataset` → display "Arm" (arm-level merge across teratomas).
  x <- ifelse(x %in% c("dataset", "Dataset"), "Arm",
              ifelse(x %in% c("celltype", "Celltype"), "Cell type",
                     stringr::str_to_title(x)))
  factor(x, levels = c("Arm", "Sample", "Cell type"))
}
set_levels <- c("Sites called\nin NTC",
                sprintf("Ground-truth sites\n(\u2265 %d%% depletion)", round(100 * TRUTH_D)))
set_title <- function(x) {
  factor(c(all_discovered = set_levels[1], pass_GT_criteria = set_levels[2])[as.character(x)],
         levels = set_levels)
}
bam_title <- function(x) {
  labs <- unname(label_bamType[as.character(x)])
  labs <- gsub(" (", "\n(", labs, fixed = T)
  factor(labs, levels = gsub(" (", "\n(", unname(label_bamType[BT_REPORTED]), fixed = T))
}

# Mismatch spectra ----
plot_bam_spectrum <- function(dt, stem, title, x_lab, fill_vals, width, height, part) {
  d <- copy(dt)
  d[, varType := factor(as.character(varType), levels = names(fill_vals))]
  d[, tier := tier_title(tier)]
  d[, set := set_title(set)]
  d[, bam_label := bam_title(bam_type)]
  d[, pct_label := paste0(round(percent*100, digits = 1), "%")]
  d[varType != "AG", pct_label := NA_character_]
  p <- ggplot(d, aes(varType, percent, fill = varType)) +
    geom_col(colour = "grey20", linewidth = 0.2, width = 0.85) +
    geom_text(aes(label = pct_label), vjust = -1, hjust = 0.1) +
    ggh4x::facet_nested(tier ~ bam_label + set) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_y_continuous(labels = function(x) round(100 * x, 1),
                       expand = expansion(mult = c(0, 0.25)),
                       breaks = seq(0, 1, 0.25),
                       limits = c(0, 1)) +
    labs(title = title, x = x_lab, y = "Percentage of sites (%)") +
    theme_manuscript(facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          panel.spacing.x = unit(0.6, "lines"),
          legend.title.position = "none")
  save_figure_sd(p, stem, width = width, height = height, parts = part, grid = "y", script = SCRIPT, src_stem = STEM)
  invisible(p)
}
plot_bam_spectrum(read_part(STEM, "spectrum12"), "GTDEF_bar_mismatch_spectrum_by_bamtype_vs_TRUE",
                  title = "Mismatch types by BAM processing: sites called in NTC vs ground-truth sites",
                  x_lab = "Mismatch type (RNA strand)", fill_vals = cmap_varType, width = 19, height = 7.2,
                  part = "spectrum12")

plot_bam_spectrum6 <- function(dt, stem, title, width, height) {
  d <- copy(dt)
  d[, varType := factor(as.character(varType), levels = varType6_levels)]
  d[, tier := tier_title(tier)]
  d[, set := set_title(set)]
  d[, bam_label := bam_title(bam_type)]
  d[, pct_label := paste0(round(percent*100, digits = 1), "%")]
  d[as.character(varType) != "AG", pct_label := NA_character_]
  p <- ggplot(d, aes(varType, percent, fill = varType)) +
    geom_col(colour = "grey20", linewidth = 0.2, width = 0.85) +
    geom_text(aes(label = pct_label), vjust = -1) +
    ggh4x::facet_nested(tier ~ bam_label + set) +
    scale_fill_manual(values = cmap_varType, guide = "none") +
    scale_y_continuous(labels = function(x) round(100 * x, 1),
                       expand = expansion(mult = c(0, 0.25)),
                       breaks = seq(0, 1, 0.25),
                       limits = c(0, 1)) +
    labs(title = title,
         x = "Mismatch type (strand complements collapsed)",
         y = "Percentage of sites (%)") +
    theme_manuscript(facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(size = 8),
          panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          panel.spacing.x = unit(0.6, "lines"))
  save_figure_sd(p, stem, width = width, height = height, parts = "spectrum6", grid = "y", script = SCRIPT,
                 src_stem = STEM)
}
plot_bam_spectrum6(read_part(STEM, "spectrum6"), "GTDEF_bar_mismatch_spectrum6_by_bamtype_vs_TRUE",
                   title = "Collapsed mismatch types by BAM processing: sites called in NTC vs ground-truth sites",
                   width = 19, height = 7.2)

# Euler: coordinate vs UMI de-duplication, one panel per tier ----
eu <- read_part(STEM, "euler")
for (tier_i in intersect(TIERS, eu$tier)) {
  combo <- setNames(eu[tier == tier_i, n_sites], eu[tier == tier_i, region])
  fit <- euler(combo)
  cols <- c(Coord = unname(cmap_proc[["De-duplicated (coordinate)"]]),
            UMI = unname(cmap_proc[["De-duplicated (UMI)"]]))
  p_eu <- plot(fit,
               fills = list(fill = cols, alpha = 0.45),
               edges = list(col = cols, lwd = 1.5),
               labels = list(fontsize = 11,
                             labels = c(Coord = "Coord dedup",
                                        UMI = "UMI dedup")),
               quantities = list(fontsize = 9),
               main = list(
                 label = paste0("Ground-truth sites called in NTC — ",
                                as.character(tier_title(tier_i)[1])),
                 fontsize = 12, font = 2))
  stem <- paste0("GTDEF_euler_TRUE_sites_by_bamtype_", tier_file_stem(tier_i))
  svg_file <- file.path(FIG_DIR, paste0(stem, ".svg"))
  grDevices::svg(svg_file, width = 6.5, height = 5.5)
  print(p_eu)
  dev.off()
  srcs <- sd_sources(SRC_DIR, STEM, "euler")
  fwrite(data.table(plot_file = basename(svg_file), source_file = repo_relative(srcs), source_exists = file.exists(srcs)),
         file.path(SRC_DIR, paste0(stem, ".sources.tsv")), sep = "\t")
}

# UpSet across the reported BAM methods (ggplot; UpSetR not in env) ----
combos <- read_part(STEM, "upset_combos")
sets <- read_part(STEM, "upset_sets")
upset_from_counts <- function(tier_i, title, max_combos = 40L) {
  comb_n <- combos[tier == tier_i, .(comb, n)]
  bt_levels <- strsplit(combos[tier == tier_i, bam_types][1], ",")[[1]]
  setorder(comb_n, -n)
  if (nrow(comb_n) > max_combos) comb_n <- comb_n[seq_len(max_combos)]
  comb_n[, idx := .I]
  decode <- function(key) as.integer(strsplit(key, "")[[1]])
  bar_dt <- comb_n[, .(idx, n, comb)]
  mat_long <- rbindlist(lapply(seq_len(nrow(comb_n)), function(i) {
    bits <- decode(comb_n$comb[i])
    data.table(idx = i, bam_type = bt_levels,
               present = factor(bits == 1L, levels = c(FALSE, TRUE)))
  }))
  mat_long[, bam_label := factor(unname(label_bamType[bam_type]),
                                 levels = unname(label_bamType[bt_levels]))]
  set_n <- sets[tier == tier_i, .(bam_type, n)]
  set_n[, bam_label := factor(unname(label_bamType[bam_type]),
                              levels = unname(label_bamType[bt_levels]))]

  p_bar <- ggplot(bar_dt, aes(idx, n)) +
    geom_col(fill = "#e6e9f5", width = 0.7, color = "black", linewidth = 0.25) +
    scale_x_continuous(breaks = bar_dt$idx, expand = expansion(add = 0.3)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       n.breaks = 6,
                       labels = shorten_axis) +
    labs(x = NULL, y = "Intersection size", title = title) +
    theme_manuscript() + grid_y +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 14, hjust = 0.5))

  p_mat <- ggplot() +
    geom_point(data = mat_long[present == FALSE],
               aes(idx, forcats::fct_rev(bam_label),
                   colour = present, size = present)) +
    geom_line(data = mat_long[present == TRUE],
              aes(idx, forcats::fct_rev(bam_label),
                  group = idx),
              colour = "grey20", linewidth = 0.4) +
    geom_point(data = mat_long[present == TRUE],
               aes(idx, forcats::fct_rev(bam_label),
                   colour = present, size = present)) +
    scale_colour_manual(values = c("FALSE" = "grey85", "TRUE" = "black"),
                        guide = "none") +
    scale_size_manual(values = c("FALSE" = 1.6, "TRUE" = 2.8), guide = "none") +
    scale_x_continuous(breaks = bar_dt$idx, expand = expansion(add = 0.6)) +
    labs(x = "Intersection", y = NULL) +
    theme_manuscript() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          panel.grid.major.x = element_blank())

  p_set <- ggplot(set_n, aes(n, forcats::fct_rev(bam_label), fill = bam_label)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.25, show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = shorten_axis) +
    scale_fill_manual(values = cmap_proc) +
    labs(x = "Set size", y = NULL) +
    theme_manuscript() + grid_x +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  egg::ggarrange(p_bar, ggplot() + theme_void(),
                 p_mat, p_set, byrow = TRUE, nrow = 2,
                 heights = c(0.6, 0.4), widths = c(0.75, 0.25))
}

UPSET_PARTS <- c("upset_combos", "upset_sets")
tiers_up <- intersect(TIERS, combos$tier)
upset_plots <- list()
for (tier_i in tiers_up) {
  wide_n <- combos[tier == tier_i, sum(n)]
  upset_plots[[tier_i]] <- upset_from_counts(
    tier_i,
    title = paste0("Ground-truth sites called in NTC by BAM processing method: ",
                   as.character(tier_title(tier_i)[1]),
                   " (n = ", format(wide_n, big.mark = ","), ")"))
}
for (tier_i in names(upset_plots)) {
  stem <- paste0("GTDEF_upset_TRUE_sites_by_bamtype_", tier_file_stem(tier_i))
  save_figure_sd(upset_plots[[tier_i]], stem, width = 10, height = 6, parts = UPSET_PARTS, grid = "y",
                 script = SCRIPT, src_stem = STEM)
}

# Combined upset panel (tiers as rows)
upset_all_plots <- lapply(tiers_up, function(tier_i) upset_from_counts(tier_i, title = as.character(tier_title(tier_i)[1])))
if (length(upset_all_plots)) {
  p_up_panel <- cowplot::plot_grid(plotlist = upset_all_plots, ncol = 1,
                                   labels = "auto", label_size = 14)
  save_figure_sd(p_up_panel, "GTDEF_upset_TRUE_sites_by_bamtype", width = 9.5, height = 12.5, parts = UPSET_PARTS,
                 grid = "y", script = SCRIPT, src_stem = STEM)
}
