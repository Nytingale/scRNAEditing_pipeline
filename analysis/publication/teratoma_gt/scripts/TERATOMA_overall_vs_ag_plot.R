#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_overall_vs_ag_plot.R   (HANDOFF 2026-09-29 UPDATE-27)
# A-to-G percentage (truth-free) as a proxy for performance: cross-tier summaries of each reported
# configuration (full benchmark call set) against its A-to-G percentage, one panel per truth set.
#   mean F1      = mean of the three tier F1 values (tier F1 from unit-averaged precision and recall)
#   overall rank = rank of the summed per-tier F1 ranks, ties broken by the higher mean F1 (as in the grids)
#   A-to-G %     = per tier, A-to-G sites / all sites pooled over the tier's units; then the mean of the
#                  three tiers (the same value in every panel)
#
# Inputs
#   the nine truth sets in TERATOMA_common.R (NINE_TRUTHS; 2026-10-01)
#   human/pm/gt_eval_v6/scoring/sweep_infneg/teratoma_inputs/ag_percent_units.tsv
# Outputs (figures/pm_gt/)
#   TERATOMA_scatter_ag_vs_mean_f1_by_truth.svg
#   TERATOMA_scatter_ag_vs_overall_rank_by_truth.svg
#   TERATOMA_scatter_ag_sites_vs_mean_f1_by_truth.svg    A-to-G sites per unit: per tier mean over units,
#                                                        geometric mean of the three tiers
#   TERATOMA_scatter_ag_composite_rank_vs_mean_f1_by_truth.svg  A-to-G composite rank = mean of the A-to-G %
#                                                        rank and the A-to-G sites rank (1 = best on each; no truth used)
# Single-truth (ground truth 50%) version: ../_superseded/overall_vs_ag_single_2026-10-01/
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_overall_vs_ag_plot.R"

COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
TRUTH_LAB <- setNames(NINE_TRUTHS$label, NINE_TRUTHS$truth)
SRC <- c(setNames(NINE_TRUTHS$file, NINE_TRUTHS$truth), ag = file.path(SCORING_DIR, "ag_percent_units.tsv"))

overall <- function(s) {
  s <- copy(s)[, f1 := fifelse(precision + recall_vs_union > 0,
                               2 * precision * recall_vs_union / (precision + recall_vs_union), NA_real_)]
  s[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
  r <- s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
  r[, overall_rank := frank(list(rank_sum, -mean_f1), ties.method = "min")]
}
ov <- rbindlist(lapply(names(TRUTH_LAB), function(n) overall(read_truth_scoring(n))[, truth := n]))

ag <- fread(SRC[["ag"]])
ag[, (COMBO) := tstrsplit(config, "__", fixed = TRUE)]
ag_ov <- ag[, .(ag_pct = 100 * sum(n_ag) / sum(n_sites), ag_sites = mean(n_ag)), by = c("tier", COMBO)][
  , .(ag_percent = mean(ag_pct), ag_sites = exp(mean(log(ag_sites + 1))) - 1, n_tiers = .N), by = COMBO]
# A-to-G composite rank: mean of the rank on A-to-G % and the rank on A-to-G sites (1 = best on each); truth-free
ag_ov[, proxy_rank := (frank(-ag_percent, ties.method = "average") + frank(-ag_sites, ties.method = "average")) / 2]

d <- merge(ov, ag_ov, by = COMBO)
stopifnot(nrow(d) == nrow(ov), all(d$n_tiers == 3L))
strip_lab <- function(x) {
  x <- gsub("depletion", "editing reduction in KO", x, fixed = TRUE)
  x <- gsub(", KO cov", ", NTC and KO cov", x, fixed = TRUE)
  x <- gsub(">= 1 G", ">= 1 variant read", x, fixed = TRUE)
  x <- gsub("([0-9]+) G\\b", "\\1 variant reads", x)
  x <- gsub("edit >=", "editing level >=", x, fixed = TRUE)
  gsub(", ", ",\n", gsub(">=", "\u2265", x, fixed = TRUE), fixed = TRUE)
}
TRUTH_STRIP <- setNames(strip_lab(TRUTH_LAB), names(TRUTH_LAB))
d[, `:=`(Truth = factor(TRUTH_STRIP[truth], TRUTH_STRIP), method = relabel(caller, label_caller),
         proc = relabel(bam_type, label_bamType), is_rank1 = overall_rank == 1L)]
PROC_LEVELS <- intersect(label_bamType[c("sampleRmDup", "sampleUmiDup", "recal", "recalq20", "recaldefault")], d$proc)
d[, proc := factor(proc, PROC_LEVELS)]
shape_proc <- setNames(c(21, 24, 22, 23, 25)[seq_along(PROC_LEVELS)], PROC_LEVELS)

# Pearson on log10 counts for the A-to-G sites axis (plotted on a log scale)
stats <- function(xvar, yvar) d[, {x <- if (xvar == "ag_sites") log10(get(xvar) + 1) else get(xvar)
  .(pearson = cor(x, get(yvar)), spearman = cor(x, get(yvar), method = "spearman"))}, by = Truth][
  , lab := sprintf("Pearson's r = %.2f\nSpearman's rho = %.2f", pearson, spearman)]

X_SPEC <- list(
  ag_percent = list(lab = "A-to-G percentage (%)", scale = scale_x_continuous(breaks = seq(20, 100, 20)), right = FALSE),
  ag_sites   = list(lab = "A-to-G sites per unit", right = FALSE,
                    scale = scale_x_log10(labels = shorten_axis, breaks = c(10, 30, 100, 300, 1000, 3000, 10000, 30000))),
  proxy_rank = list(lab = "A-to-G composite rank", right = TRUE,
                    scale = scale_x_continuous(breaks = c(1, seq(10, 50, 10)))))

# shape_by_proc: shape = preprocessing (filled shapes in the caller colour);
# otherwise shape marks overall rank 1 (triangle)
scatter <- function(yvar, ylab, title, stem, xvar = "ag_percent", shape_by_proc = FALSE) {
  st <- stats(xvar, yvar)
  up <- yvar == "mean_f1"
  xs <- X_SPEC[[xvar]]
  lab_right <- xs$right == up
  pts <- if (shape_by_proc) list(
    geom_point(aes(fill = method, shape = proc), colour = "grey15", stroke = 0.3, size = 4.2, alpha = 0.6),
    scale_fill_manual(values = cmap_method, name = NULL,
                      guide = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, size = 4.5, alpha = 0.9))),
    scale_shape_manual(values = shape_proc, name = NULL,
                       guide = guide_legend(nrow = 1, order = 2, override.aes = list(fill = "grey60", size = 4.5, alpha = 0.9)))
  ) else list(
    geom_point(aes(colour = method, shape = is_rank1), size = 2.6, alpha = 0.85),
    scale_colour_manual_cvd(values = cmap_method, name = NULL,
                            guide = guide_legend(nrow = 1, override.aes = list(shape = 16, size = 3.5))),
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none"))
  if (shape_by_proc) assert_cvd_safe(cmap_method, "scale_fill_manual")
  p <- ggplot(d, aes(.data[[xvar]], .data[[yvar]])) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    pts +
    geom_text(data = st, aes(x = if (lab_right) max(d[[xvar]]) else min(d[[xvar]]), y = Inf, label = lab),
              hjust = if (lab_right) 1 else 0, vjust = if (up) 1.3 else -0.3, size = 4.6, lineheight = 1.05) +
    facet_wrap(~ Truth, nrow = 3, scales = if (up) "free_y" else "fixed") +
    xs$scale +
    labs(x = xs$lab, y = ylab, title = title) +
    theme_manuscript(base_size = 13, facets = TRUE) + grid_both + theme(aspect.ratio = 1)
  if (shape_by_proc) p <- p + theme(legend.box = "vertical", legend.spacing.y = unit(2, "pt"))
  p <- if (up) p + scale_y_continuous(breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.05, 0.22)))
       else p + scale_y_reverse(breaks = c(1, seq(10, 50, 10)), expand = expansion(mult = c(0.05, 0.05)))
  save_figure(p, stem, width = 12, height = 14, sources = SRC, data = d, grid = "both", script = SCRIPT)
  st[, .(Truth, pearson = round(pearson, 3), spearman = round(spearman, 3))]
}
print(scatter("mean_f1", "Mean F1 across tiers", "Mean F1 vs A-to-G percentage", "TERATOMA_scatter_ag_vs_mean_f1_by_truth",
              shape_by_proc = TRUE))
within <- d[, .(n = .N, spearman = cor(ag_percent, mean_f1, method = "spearman")), by = .(Truth, caller)]
cat("\nSpearman (A-to-G vs mean F1) within caller:\n"); print(dcast(within, caller ~ Truth, value.var = "spearman")[], digits = 2)
print(scatter("overall_rank", "Overall cross-tier rank", "Overall rank vs A-to-G percentage", "TERATOMA_scatter_ag_vs_overall_rank_by_truth"))
print(scatter("mean_f1", "Mean F1 across tiers", "Mean F1 vs number of A-to-G sites",
              "TERATOMA_scatter_ag_sites_vs_mean_f1_by_truth", xvar = "ag_sites", shape_by_proc = TRUE))
print(scatter("mean_f1", "Mean F1 across tiers", "Mean F1 vs A-to-G composite rank",
              "TERATOMA_scatter_ag_composite_rank_vs_mean_f1_by_truth", xvar = "proxy_rank", shape_by_proc = TRUE))
