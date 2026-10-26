#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_results_plot.R        (stage 2 of 2 -- FIGURES)
#
# The teratoma ground-truth benchmark RESULT. The panels that justify the
# ground-truth definition are a separate deliverable (publication/gt_definition).
#
# >>> THE UNIT OF COMPARISON IS THE UNIQUE METHOD COMBINATION. <<<
# caller x bam_type x strand_mode x snp_subtract, never caller alone: "how does
# reditools3 do" has no answer without stating the dedup, the stranding and the
# SNP treatment. ALL 64 combinations are shown, faceted by caller so each panel
# sizes to the combinations that caller actually has.
#
# Axis text uses the manuscript's descriptive names (label_bamType etc.), never
# raw pipeline tokens: a reader has no way to know what "sampleRmDup" is.
#
# Panels
#   (a) grouped bars  precision / recall / F1, all 64 combinations
#   (b) scatter       discovered sites vs F1 -- best is upper right
#   (c) scatter       discovered sites vs A-to-G percentage (ground-truth-free)
#   (d) combined figure assembling the three
#
# Usage: Rscript TERATOMA_results_plot.R [DATASET] [TIER]
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "TERATOMA_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
teratoma_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
TIER   <- if (length(argv) >= 2) argv[2] else "sample"
SCRIPT <- caller_script()
SRC    <- file.path(REPORT_DIR, "TERATOMA_method_performance.tsv")
if (!file.exists(SRC)) stop("missing ", SRC, " -- run TERATOMA_prep_source_data.R first")

d <- fread(SRC)[tier == TIER]
if (!nrow(d)) stop("no rows for tier '", TIER, "'")

# ── Per-replicate scores for jitter points ──────────────────────────────────
# sample tier: one point per teratoma (T1–T4).
# celltype tier: one point per cell type (17 cell types).
# dataset tier: no biological replicates — crossbar only, no jitter.
REP_TSV <- file.path(SCORING_DIR, "gt_tier_replicate_scores.tsv")
HAVE_REPS <- file.exists(REP_TSV) && TIER %in% c("sample", "celltype")
if (HAVE_REPS) {
  reps <- fread(REP_TSV)[tier == TIER]
  reps[, f1 := fifelse(precision + recall_vs_union > 0,
                       2 * precision * recall_vs_union / (precision + recall_vs_union),
                       NA_real_)]
  if (TIER == "sample") {
    reps[, replicate := teratoma]
    REP_CMAP  <- cmap_teratoma_jitter
    REP_LABEL <- "Teratoma"
  } else {
    reps[, replicate := celltype]
    REP_CMAP  <- cmap_teratoma_celltype
    REP_LABEL <- "Cell type"
  }
  message("  replicate scores: ", nrow(reps), " rows, ",
          uniqueN(reps$replicate), " ", tolower(REP_LABEL), "s")
} else {
  reps <- NULL
  REP_CMAP <- NULL; REP_LABEL <- NULL
}

# Official method spellings on every plot, and the manuscript's descriptive
# grid wording, never raw pipeline tokens.
# >>> FIGURE STEMS INCLUDE THE TIER. <<<
# Running sample then dataset used to overwrite the same SVG names. Every
# stem below is suffixed with _${tier_file_stem(TIER)} (`dataset` → `arm`).
TIER_STEM <- tier_file_stem(TIER)
d[, method := relabel(caller, label_caller)]
# TWO strings, deliberately. The DISPLAYED label omits the method, because
# colour already encodes it -- but bam_type | strand | snp is not unique across
# methods, so using it as the factor level would silently MERGE rows that share
# it. The factor is keyed on the full combination and only the display text is
# shortened.
d[, lab_key := paste(method,
                     relabel(bam_type, label_bamType),
                     relabel(strand_mode, label_strandMode),
                     relabel(snp_subtract, label_snpSubtract), sep = " | ")]
d[, lab_show := paste(relabel(bam_type, label_bamType),
                      relabel(strand_mode, label_strandMode),
                      relabel(snp_subtract, label_snpSubtract), sep = "\n")]
stopifnot(!anyDuplicated(d$lab_key))
setorder(d, -f1)

# Points overlap heavily where combinations score alike; without alpha the
# clusters read as single points and hide how tightly the leaders are packed.
PT_ALPHA <- 0.65
PT_SIZE  <- 2.8
# Scatters: higher alpha so overlapping methods stay readable; rank-1 (best F1
# in this tier) is drawn as a triangle on the A-to-G panels.
SCATTER_ALPHA <- 0.7
SCATTER_SIZE  <- 2.8
# Scatters are square: both axes are "more is better", so a 1:1 box keeps the
# distance from the upper-right corner visually honest.
sq <- theme(aspect.ratio = 1)

scatter_ag_base <- function(dat, xvar, xlab, title) {
  dd <- copy(dat)
  dd[, is_rank1 := frank(-f1, ties.method = "min") == 1L]
  dd[, xval := get(xvar)]
  ggplot(dd, aes(xval, ag_percent, colour = method, shape = is_rank1)) +
    geom_hline(yintercept = 100 / 12, linetype = 2, colour = "black",
               linewidth = 0.4) +
    annotate("text", x = max(dd$xval, na.rm = TRUE), y = 100 / 12,
             hjust = 1, vjust = 1.6, size = 3.4, colour = "black",
             label = "Chance level if mismatches were random (1 in 12)") +
    geom_point(size = SCATTER_SIZE, alpha = SCATTER_ALPHA) +
    scale_colour_manual_cvd(values = cmap_method, name = NULL,
                            guide = guide_legend(nrow = 2, override.aes = list(shape = 16))) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = xlab, y = "A-to-G percentage (%)", title = title) +
    theme_manuscript() + grid_both + sq
}

# (a) precision / recall / F1, ALL combinations ----
# FACET BY METRIC, COLOUR BY METHOD. One row per combination appearing once in
# each of three metric columns, rather than three dodged bars per row: at 64
# combinations the dodged form is 192 bars and the eye cannot follow a single
# combination across metrics. Colour now carries the method, which is the
# comparison the panel is for.
STEM <- paste0("TERATOMA_bar_method_performance_", TIER_STEM)
idv <- c("lab_key", "lab_show", "method", "f1")
dl <- melt(d, id.vars = idv,
           measure.vars = c("precision", "recall_vs_union", "f1"),
           variable.name = "metric", value.name = "value")
dl[, metric := factor(metric, levels = c("precision", "recall_vs_union", "f1"),
                      labels = c("Precision", "Recall", "F1"))]
# facet_grid(preprocessing ~ metric), one bar per method inside each cell.
#
# The shared bam_type | strand | snp label is the GROUPING, not a collision:
# putting the methods together in one cell is what makes them comparable at a
# fixed preprocessing configuration, which is the comparison the panel is for.
#
# The method names are NOT repeated on the discrete axis -- colour and the
# legend already identify them, and 5 labels x 24 rows is pure noise. Strip
# labels are three lines rather than pipe-separated, since a y-side strip has
# the width for it.
ALL_METHODS <- sort(unique(d$method))

build_jitter_panel <- function(dat, title, rep_data = NULL,
                               rep_cmap = NULL, rep_label = NULL,
                               show_rank_markers = TRUE) {
  # Complete the method × config grid so capability gaps show "NA" rather than
  # a score of zero (see main.nf: hybrid stranding, snp_subtract).
  full <- CJ(lab_show = unique(dat$lab_show), metric = unique(dat$metric),
             method = ALL_METHODS, unique = TRUE)
  full <- merge(full, dat[, .(lab_show, metric, method, value)],
                by = c("lab_show", "metric", "method"), all.x = TRUE)
  full[, missing := !is.finite(value)]

  # Replicate jitter data (per-teratoma or per-celltype).
  have_jitter <- !is.null(rep_data) && nrow(rep_data) > 0
  if (have_jitter) {
    jd <- melt(rep_data,
               id.vars = c("lab_key", "lab_show", "method", "replicate"),
               measure.vars = c("precision", "recall_vs_union", "f1"),
               variable.name = "metric", value.name = "value")
    jd[, metric := factor(metric, levels = c("precision", "recall_vs_union", "f1"),
                          labels = c("Precision", "Recall", "F1"))]
  }

  # When replicate data exists, use the replicate MEAN for row ordering,
  # top-3 ranking and crossbar position — a method is best when it is
  # reproducibly and robustly best, not when its pooled-union score is high.
  if (have_jitter) {
    rep_means <- jd[is.finite(value),
                    .(value = mean(value, na.rm = TRUE)),
                    by = .(lab_show, metric, method)]
  } else {
    rep_means <- NULL
  }
  score_for_rank <- if (!is.null(rep_means)) rep_means else dat

  # Best F1 first (top row).
  ord <- score_for_rank[metric == "F1",
           .(best = max(value, na.rm = TRUE)),
           by = lab_show][order(-best), as.character(lab_show)]
  full[, lab_show := factor(as.character(lab_show), levels = ord)]
  full[, method := factor(method, levels = rev(ALL_METHODS))]
  
  if (have_jitter) {
    jd[, lab_show := factor(as.character(lab_show), levels = ord)]
    jd[, method := factor(method, levels = rev(ALL_METHODS))]
    rep_means[, method := factor(method, levels = levels(full$method))]
    jd[, replicate := factor(replicate, levels = names(rep_cmap))]
  }

  # Top-3 F1 markers — ranked by replicate mean when available.
  rank_src <- if (!is.null(rep_means)) rep_means else full
  top3 <- rank_src[metric == "F1" & is.finite(value)][order(-value)][seq_len(min(3, .N))]
  top3[, rank := seq_len(.N)]
  # Place marker past the highest replicate point, not the mean.
  if (have_jitter) {
    rep_max <- jd[metric == "F1" & is.finite(value),
                  .(max_val = max(value, na.rm = TRUE)),
                  by = .(lab_show, method)]
    top3 <- merge(top3, rep_max, by = c("lab_show", "method"), all.x = TRUE)
    top3[, mark_x := pmin(1, fifelse(is.finite(max_val), max_val, value)) + 0.035]
    top3[, max_val := NULL]
  } else {
    top3[, mark_x := pmin(1, value) + 0.04]
  }
  top3[, lab_show := factor(as.character(lab_show), levels = ord)]

  p <- ggplot(full, aes(value, method))

  # Jitter layer: one point per replicate, coloured by teratoma or cell type.
  # Only y-jitter (method axis); x stays at the exact score position.
  # Points (no jitter, exact position) then transparent boxplots on top.
  if (have_jitter) {
    p <- p +
      geom_point(data = jd, aes(colour = replicate),
                 size = 1.4, alpha = 0.5, na.rm = TRUE) +
      geom_boxplot(data = jd, aes(group = method),
                   fill = NA, colour = "black", linewidth = 0.35,
                   width = 0.45, outlier.shape = NA, na.rm = TRUE) +
      scale_colour_manual(values = rep_cmap, name = rep_label,
                          guide = guide_legend(
                            title.position = "top", ncol = 3,
                            override.aes = list(size = 2.5, alpha = 1)))
  } else {
    # Dataset tier: no replicates — show bars coloured by method.
    p <- p +
      geom_col(data = full[is.finite(value)], aes(fill = method),
               orientation = "y", width = 0.72, colour = "grey20",
               linewidth = 0.15, na.rm = TRUE) +
      scale_fill_manual_cvd(values = cmap_method, name = NULL, guide = "none")
  }

  p <- p +
    geom_text(data = full[missing == TRUE], aes(x = 0.015, y = method),
              label = "NA", hjust = 0, size = 2.3, colour = "black",
              inherit.aes = FALSE)

  if (show_rank_markers) {
    p <- p +
      geom_text(data = top3, aes(x = mark_x, y = method),
                label = "\u25C0", hjust = 0, size = 3.1, colour = "black",
                inherit.aes = FALSE) +
      geom_text(data = top3, aes(x = mark_x + 0.07, y = method, label = rank),
                hjust = 0, size = 3.2, fontface = "bold", colour = "black",
                inherit.aes = FALSE)
  }

  p <- p +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.02, 1.08),
                       expand = expansion(mult = c(0, 0))) +
    facet_grid(lab_show ~ metric, switch = "y") +
    labs(x = "Score", y = NULL, title = title) +
    theme_manuscript(facets = TRUE) + grid_x + strip_y_outside +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9),
          panel.spacing.x = unit(1.5, "lines"),
          panel.spacing.y = unit(0.55, "lines"),
          axis.text.y = element_text(size = 8),
          legend.position = if (have_jitter) "bottom" else "none",
          legend.justification = "left",
          legend.box.just = "left",
          legend.key.height = unit(0.7, "lines"),
          legend.spacing.y = unit(0.15, "lines"))
  p
}

# >>> SPLIT ON HOW MANY METHODS COULD RUN THE CONFIGURATION. <<<
# A facet is a CROSS-METHOD comparison, so one with fewer than three methods
# is mostly NA and compares almost nothing. Those go to a supplement.
#
# The cost is quantified rather than assumed: among the 64 bars, those in
# <3-method facets top out at F1 0.575 (BCFtools | Unprocessed | Hybrid |
# + blacklist), which is rank 5 of 64, behind four annotation bars at 0.592,
# 0.591, 0.588 and 0.581. So the supplement holds nothing in the top four, and
# the main figure keeps every bar that is genuinely in contention.
#
# In practice this moves the hybrid-stranding facets (BCFtools + GATK only,
# since hybrid needs I16/SB) and the recaldefault facets (GATK only).
MIN_METHODS_MAIN <- 3L
dl[, n_methods_in_facet := uniqueN(method[is.finite(value)]), by = lab_show]
main <- dl[n_methods_in_facet >= MIN_METHODS_MAIN]
supp <- dl[n_methods_in_facet <  MIN_METHODS_MAIN]

TIER_LABEL <- c(sample = "Sample level", dataset = "Arm level",
                celltype = "Cell type level")
tier_lab <- TIER_LABEL[TIER]
if (is.na(tier_lab)) tier_lab <- paste0(TIER, " tier")

n_rows <- uniqueN(main$lab_show)
# Prepare replicate data with the same label keys as the pooled data.
main_rep <- NULL
if (HAVE_REPS && !is.null(reps)) {
  reps[, method := relabel(caller, label_caller)]
  reps[, lab_key := paste(method,
                          relabel(bam_type, label_bamType),
                          relabel(strand_mode, label_strandMode),
                          relabel(snp_subtract, label_snpSubtract), sep = " | ")]
  reps[, lab_show := paste(relabel(bam_type, label_bamType),
                           relabel(strand_mode, label_strandMode),
                           relabel(snp_subtract, label_snpSubtract), sep = "\n")]
  main_rep <- reps[lab_show %in% unique(main$lab_show)]
}

p.a <- build_jitter_panel(
  main,
  paste0("Method comparisons at each processing configuration: ", tier_lab),
  rep_data = main_rep, rep_cmap = REP_CMAP, rep_label = REP_LABEL)
save_figure(p.a, STEM, width = 9.5, height = 2.2 + 0.75 * n_rows,
            sources = SRC, data = main, grid = "x", script = SCRIPT)

if (nrow(supp)) {
  supp_rep <- if (!is.null(main_rep)) reps[lab_show %in% unique(supp$lab_show)] else NULL
  p.supp <- build_jitter_panel(
    supp,
    paste0("Configurations fewer than ", MIN_METHODS_MAIN,
           " methods can run: ", tier_lab),
    rep_data = supp_rep, rep_cmap = REP_CMAP, rep_label = REP_LABEL,
    show_rank_markers = FALSE)
  save_figure(p.supp, paste0("TERATOMA_bar_few_methods_", TIER_STEM), width = 9.5,
              height = 2.2 + 0.75 * uniqueN(supp$lab_show),
              sources = SRC, data = supp, grid = "x", script = SCRIPT)
}

# (b) discovered sites vs F1 ----
STEM <- paste0("TERATOMA_scatter_discovery_vs_f1_", TIER_STEM)
p.b <- ggplot(d, aes(n_discovered, f1, colour = method)) +
  geom_point(size = PT_SIZE, alpha = PT_ALPHA) +
  scale_colour_manual_cvd(values = cmap_method, name = NULL,
                          guide = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = shorten_axis) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(x = "Discovered A-to-G sites", y = "F1",
       title = paste0("Yield against quality: ", tier_lab)) +
  theme_manuscript() + grid_both + sq
save_figure(p.b, STEM, width = 5.4, height = 5.3, sources = SRC, data = d,
            grid = "both", script = SCRIPT)

# (c) discovered sites vs A-to-G percentage ----
# A-to-G percentage is an INDEPENDENT specificity proxy: it needs no ground
# truth at all, so agreement with the F1 ranking is corroboration from a
# direction the truth set cannot influence.
# (d) F1 vs A-to-G percentage — same proxy against the primary quality axis.
p.c <- NULL
p.d <- NULL
if ("ag_percent" %in% names(d) && any(is.finite(d$ag_percent))) {
  da <- d[is.finite(ag_percent)]
  STEM <- paste0("TERATOMA_scatter_discovery_vs_ag_", TIER_STEM)
  p.c <- scatter_ag_base(da, "n_discovered", "Discovered A-to-G sites",
                         paste0("A-to-G enrichment vs yield: ", tier_lab)) +
    scale_x_continuous(labels = shorten_axis)
  save_figure(p.c, STEM, width = 5.4, height = 5.5, sources = SRC, data = da,
              grid = "both", script = SCRIPT)

  STEM <- paste0("TERATOMA_scatter_f1_vs_ag_", TIER_STEM)
  # Each axis runs from 0 to just above its own maximum, so the spread of the
  # points reflects the correlation rather than the unused range of a fixed
  # 0-1 / 0-100 box. The least-squares line is the fit behind Pearson's r.
  da2 <- copy(da)
  da2[, is_rank1 := frank(-f1, ties.method = "min") == 1L]
  r_pears <- cor(da2$f1, da2$ag_percent, use = "complete.obs", method = "pearson")
  r_spear <- cor(da2$f1, da2$ag_percent, use = "complete.obs", method = "spearman")
  x_max <- max(da2$f1, na.rm = TRUE) * 1.05
  y_max <- max(da2$ag_percent, na.rm = TRUE) * 1.05
  p.d <- ggplot(da2, aes(f1, ag_percent, colour = method, shape = is_rank1)) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
                linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    geom_point(size = SCATTER_SIZE, alpha = SCATTER_ALPHA) +
    scale_colour_manual_cvd(values = cmap_method, name = NULL,
                            guide = guide_legend(nrow = 2, override.aes = list(shape = 16))) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_x_continuous(breaks = scales::breaks_pretty(5), expand = expansion(mult = c(0.02, 0))) +
    scale_y_continuous(breaks = scales::breaks_pretty(5), expand = expansion(mult = c(0.02, 0))) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    annotate("text", x = Inf, y = -Inf,
             label = sprintf("Pearson's r = %.3f\nSpearman's \u03c1 = %.3f",
                             r_pears, r_spear),
             hjust = 1.08, vjust = -0.4, size = 3.6, lineheight = 1.05) +
    labs(x = "F1", y = "A-to-G percentage (%)",
         title = paste0("A-to-G enrichment vs F1: ", tier_lab)) +
    theme_manuscript() + grid_both +
    theme(aspect.ratio = 1)
  save_figure(p.d, STEM, width = 5.4, height = 5.5, sources = SRC, data = da2,
              grid = "both", script = SCRIPT)
}

# combined figure ----
# The all-combinations bar chart leads as panel (a): it is the result, and the
# two scatters are summaries of it. The square scatters share the bottom row as
# (b) and (c). rel_heights comes from the bar panel's row count, so adding a
# combination lengthens the figure rather than squashing the chart.
if (!is.null(p.c)) {
  bar_h <- 2.2 + 0.75 * n_rows
  # Panel (a) uses a narrower width so the y-axis strip labels (preprocessing
  # configurations) are not clipped. The scatters share the bottom row.
  row_bar <- cowplot::plot_grid(p.a, ncol = 1, labels = "a", label_size = 18)
  row_sc  <- cowplot::plot_grid(p.b, p.c, ncol = 2, labels = c("b", "c"),
                                label_size = 18)
  p.panel <- cowplot::plot_grid(row_bar, row_sc, ncol = 1,
                                rel_heights = c(bar_h, 5.6))
  save_figure(p.panel, paste0("TERATOMA_panel_results_", TIER_STEM), width = 9.5,
              height = bar_h + 5.6, sources = c(SRC), data = NULL,
              grid = "both", rasterized = FALSE, script = SCRIPT)
}

message("\nresults figures -> ", FIG_DIR)
