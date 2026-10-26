#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_overall_vs_ag_plot.R   (2026-10-13; from TERATOMA_overall_vs_ag_plot.R)
# A-to-G percentage (truth-free) as a proxy for performance on the v8 truths and call state:
# cross-tier summaries of each reported configuration against its A-to-G percentage, one panel per
# v8 truth set (definition: human/pm/gt_eval_v8/README.md; paths in plotting/gt_v8_paths.R).
#   mean F1      = mean of the three tier F1 values
#   overall rank = rank of the summed per-tier F1 ranks, ties broken by the higher mean F1 (as in the grids)
#   A-to-G %     = per tier, A-to-G sites / all sites pooled over the tier's units, in the mismatch-state
#                  call sets (*_step5_nomm, the sets that are scored); then the mean of the three tiers
#
# Inputs
#   human/pm/gt_eval_v8/compare_config_rank.tsv (KO truths, state mismatch)
#   human/pm/gt_eval_v8/rediportal/scoring_mismatch/gt_tier_scoring_*.tsv
#   human/pm/gt_eval_v8/scoring/mismatch/ag_percent_units.tsv (human/scripts/gt_v8_ag_percent.R)
# Outputs (figures/pm_gt/)
#   TERATOMA_v8_scatter_ag_vs_mean_f1_by_truth.svg
#   TERATOMA_v8_scatter_ag_sites_vs_mean_f1_by_truth.svg    A-to-G sites per unit: per tier mean over units,
#                                                           geometric mean of the three tiers
#   TERATOMA_v8_scatter_ag_composite_rank_vs_mean_f1_by_truth.svg  A-to-G composite rank = mean of the A-to-G %
#                                                           rank and the A-to-G sites rank (1 = best on each; no truth used)
# Supersedes TERATOMA_scatter_ag_*_by_truth (v6 truths, final call sets).
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_overall_vs_ag_plot.R"

COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
TRUTH_LAB <- V8_TRUTH_LABELS[V8_TRUTH_IDS]
REDI_FILES <- v8_redi_files(V8_STATE)
SRC <- c(ko = V8_CMP_RANK, REDI_FILES, ag = file.path(v8_scoring_dir(V8_STATE), "ag_percent_units.tsv"))
stopifnot(all(file.exists(SRC)))

overall <- function(s) {
  s <- copy(s)   # f1 = mean of the unit F1 (gt_score_rediportal_v8.R), as the KO truths
  s[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
  r <- s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
  r[, overall_rank := frank(list(-mean_f1, rank_sum), ties.method = "min")]
}
ov_ko <- fread(V8_CMP_RANK)[state == V8_STATE, c(COMBO, "rank_sum", "mean_f1", "overall_rank", "gate", "t"), with = FALSE]
ov_ko[, truth := paste(gate, t, sep = "_")][, c("gate", "t") := NULL]
ov <- rbind(ov_ko, rbindlist(lapply(names(REDI_FILES), function(n) overall(fread(REDI_FILES[[n]]))[, truth := n])))
stopifnot(setequal(ov$truth, V8_TRUTH_IDS), ov[, .N, by = truth][, all(N == 58L)])

ag <- fread(SRC[["ag"]])
ag[, (COMBO) := tstrsplit(config, "__", fixed = TRUE)]
ag_ov <- ag[, .(ag_pct = 100 * sum(n_ag) / sum(n_sites), ag_sites = mean(n_ag)), by = c("tier", COMBO)][
  , .(ag_percent = mean(ag_pct), ag_sites = exp(mean(log(ag_sites + 1))) - 1, n_tiers = .N), by = COMBO]
# A-to-G composite rank: mean of the rank on A-to-G % and the rank on A-to-G sites (1 = best on each); truth-free,
# ranked among the 58 ranked configurations only
ag_ov <- ag_ov[unique(ov[, ..COMBO]), on = COMBO, nomatch = NULL]
stopifnot(nrow(ag_ov) == 58L)
ag_ov[, proxy_rank := (frank(-ag_percent, ties.method = "average") + frank(-ag_sites, ties.method = "average")) / 2]

d <- merge(ov, ag_ov, by = COMBO)
stopifnot(nrow(d) == nrow(ov), all(d$n_tiers == 3L))
strip_lab <- function(x) gsub(", ", ",\n", x, fixed = TRUE)
TRUTH_STRIP <- setNames(strip_lab(TRUTH_LAB), names(TRUTH_LAB))
d[, `:=`(Truth = factor(TRUTH_STRIP[truth], TRUTH_STRIP), method = caller_factor(caller),
         proc = relabel(bam_type, label_bamType), is_rank1 = overall_rank == 1L)]
PROC_LEVELS <- intersect(label_bamType[c("sampleRmDup", "sampleUmiDup", "recal", "recalq20", "recaldefault")], d$proc)
d[, proc := factor(proc, PROC_LEVELS)]
shape_proc <- setNames(c(21, 24, 22, 23, 25)[seq_along(PROC_LEVELS)], PROC_LEVELS)

# Pearson on log10 counts for the A-to-G sites axis (plotted on a log scale)
stats <- function(xvar, yvar) d[, {x <- if (xvar == "ag_sites") log10(get(xvar) + 1) else get(xvar)
  .(pearson = cor(x, get(yvar)), spearman = cor(x, get(yvar), method = "spearman"))}, by = Truth][
  , lab := sprintf("Pearson's r = %.2f\n%s", pearson, spearman_lab(spearman))]

X_SPEC <- list(
  ag_percent = list(lab = "A-to-G percentage (%)", scale = scale_x_continuous(breaks = seq(20, 100, 20)), right = FALSE),
  ag_sites   = list(lab = "A-to-G sites per unit", right = FALSE,
                    scale = scale_x_log10(labels = shorten_axis, breaks = c(10, 30, 100, 300, 1000, 3000, 10000, 30000))),
  proxy_rank = list(lab = "A-to-G composite rank", right = TRUE,
                    scale = scale_x_continuous(breaks = c(1, seq(10, 50, 10), 58), limits = c(1, 58))))

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
    labs(x = xs$lab, y = ylab, title = title, subtitle = "Call sets: sites near other mismatch types removed") +
    theme_manuscript(base_size = 13, facets = TRUE) + grid_both + theme(aspect.ratio = 1)
  if (shape_by_proc) p <- p + theme(legend.box = "vertical", legend.spacing.y = unit(2, "pt"))
  p <- if (up) p + scale_y_continuous(breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.05, 0.22)))
       else p + scale_y_reverse(breaks = c(1, seq(10, 50, 10), 58), limits = c(58, 1), expand = expansion(mult = c(0.05, 0.05)))
  save_figure(p, stem, width = 12, height = 14, sources = SRC, data = d, grid = "both", script = SCRIPT)
  st[, .(Truth, pearson = round(pearson, 3), spearman = round(spearman, 3))]
}
print(scatter("mean_f1", "Mean F1 across tiers", "Mean F1 vs A-to-G percentage", "TERATOMA_v8_scatter_ag_vs_mean_f1_by_truth",
              shape_by_proc = TRUE))
within <- d[, .(n = .N, spearman = cor(ag_percent, mean_f1, method = "spearman")), by = .(Truth, caller)]
cat("\nSpearman (A-to-G vs mean F1) within caller:\n"); print(dcast(within, caller ~ Truth, value.var = "spearman")[], digits = 2)
print(scatter("mean_f1", "Mean F1 across tiers", "Mean F1 vs number of A-to-G sites",
              "TERATOMA_v8_scatter_ag_sites_vs_mean_f1_by_truth", xvar = "ag_sites", shape_by_proc = TRUE))
print(scatter("mean_f1", "Mean F1 across tiers", "Mean F1 vs A-to-G composite rank",
              "TERATOMA_v8_scatter_ag_composite_rank_vs_mean_f1_by_truth", xvar = "proxy_rank", shape_by_proc = TRUE))

# ---- review figure: why a composite (a: rho of mean F1 with each measure, all nine truths;
#      b-c: mean F1 vs A-to-G % and vs A-to-G sites under one truth) ----------------------------------------
EX_TRUTH <- "expected_t50"
MEAS <- c(pct = "A-to-G percentage", sites = "A-to-G sites", comp = "A-to-G composite\n(rank reversed)")
cl <- d[, .(pct = cor(ag_percent, mean_f1, method = "spearman"), sites = cor(ag_sites, mean_f1, method = "spearman"),
            comp = cor(-proxy_rank, mean_f1, method = "spearman")), by = truth]
cl <- melt(cl, id.vars = "truth", variable.name = "m", value.name = "rho")
cl[, `:=`(Measure = factor(MEAS[as.character(m)], MEAS),
          Family = factor(fifelse(startsWith(truth, "redi_"), "REDIportal", "KO ground truths"),
                          c("KO ground truths", "REDIportal")))]
print(cl[, .(min = round(min(rho), 2), max = round(max(rho), 2)), by = Measure])
p_rho <- ggplot(cl, aes(Measure, rho, fill = Family)) +
  geom_point(shape = 21, size = 3.8, colour = "grey15", alpha = 0.85,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  scale_fill_manual(values = c("KO ground truths" = "#3182BD", REDIportal = "#E6550D"), name = NULL,
                    guide = guide_legend(ncol = 1, order = 3)) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = NULL, y = "Spearman's \u03C1 with mean F1\n(one point per truth set)") +
  theme_manuscript(base_size = 13) + grid_y + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ex <- d[truth == EX_TRUTH]
one <- function(xvar, xlab, xscale) {
  rho <- cor(ex[[xvar]], ex$mean_f1, method = "spearman")
  ggplot(ex, aes(.data[[xvar]], mean_f1)) +
    geom_point(aes(fill = method, shape = proc), colour = "grey15", stroke = 0.3, size = 3.8, alpha = 0.7) +
    scale_fill_manual(values = cmap_method, name = NULL,
                      guide = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, size = 4.5, alpha = 0.9))) +
    scale_shape_manual(values = shape_proc, name = NULL,
                       guide = guide_legend(nrow = 1, order = 2, override.aes = list(fill = "grey60", size = 4.5, alpha = 0.9))) +
    xscale + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = xlab, y = "Mean F1 across tiers", title = spearman_lab(rho)) +
    theme_manuscript(base_size = 13) + grid_both + theme(aspect.ratio = 1, plot.title = element_text(size = 13))
}
p_pct <- one("ag_percent", "A-to-G percentage (%)", X_SPEC$ag_percent$scale)
p_sites <- one("ag_sites", "A-to-G sites per unit (log)", X_SPEC$ag_sites$scale)
p_why <- (p_rho | p_pct | p_sites) + plot_layout(guides = "collect", widths = c(0.8, 1, 1)) +
  plot_annotation(tag_levels = "a",
                  subtitle = paste0("b, c: ", gsub("\n", " ", TRUTH_STRIP[[EX_TRUTH]]),
                                    ". Call sets: sites near other mismatch types removed")) &
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(2, "pt"),
        plot.tag = element_text(face = "bold", size = 18))
save_figure(p_why, "TERATOMA_v8_ag_f1_correlations", width = 16, height = 7, sources = SRC, data = d,
            grid = "both", script = SCRIPT)
