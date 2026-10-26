#!/usr/bin/env Rscript

# WASP reference-bias supplement: every panel of the analysis, for one dataset (draws from source data only
# since 2026-10-25).
#
# THE ARGUMENT THE PANELS MAKE, in order:
#   (A) the correction WASP offers is negligible          -> hist of dEL at 0
#   (B) and it discards most reads at editing sites       -> % removed by allele
#   (C) so most measurable sites stop being measurable    -> site attrition
# i.e. WASP is reported as a DIAGNOSTIC that bounds reference bias, and is
# deliberately NOT part of the benchmark: the cost dwarfs the correction.
#
# Every panel is faceted by coverage floor (10 / 20 / 50). The dependence on
# the floor is itself a result: the removed-read gap collapses from 1.081 pp at
# cov >= 10 to 0.034 pp at cov >= 50, where it is no longer distinguishable
# from zero (95% CI -0.023 to +0.091, p = 0.24). That STRENGTHENS "negligible" —
# wherever editing is measured confidently, the bias is undetectable.
#
# Input : source_data/<ds>/WASP_refbias.<part>.tsv (WASP_refbias_source_data.R)
# Writes, into figures/<ds>/ (SVG only; all TSVs go to source_data/<ds>/):
#   (A) WASP_hist_delta_editing_level
#   (B) WASP_bar_removed_read_alleles
#   (C) WASP_bar_site_attrition
#   (*) WASP_panel_refbias
#
# Stems are identical across datasets and the folder says which one, so the
# same panel from two platforms stays directly comparable.
#
# Usage: Rscript WASP_refbias_plot.R [dataset]   (default jzlxv3)
#        Known datasets: scripts/WASP_datasets.tsv

source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "WASP_plot_common.R"))
script_id = "publication/wasp_refbias/scripts/WASP_refbias_plot.R"
args = commandArgs(trailingOnly = TRUE)
wasp_plot_dataset(if (length(args) >= 1) args[1] else "jzlxv3")
SD = "WASP_refbias"

# In-panel annotation size, one value for all three panels so the stacked
# figure reads as one object. Facet strip size is +1 on the axis text and is
# set by theme_manuscript(), not here.
ANNOT_SIZE = 4.0
BASE_SIZE = 14

# Top-of-facet annotation, black, clear of the bars, same top expansion in
# every panel so the three agree when stacked.
#
# TWO vjust VALUES, because vjust is in units of the TEXT BLOCK's height, not
# of the panel: past vjust = 1 the text drops by (vjust - 1) x its own height.
# So one value does not give one gap. (A)'s block is 3 lines and (B)/(C)'s are
# 1, and matching the gap at 0.75 line-heights below the panel top means
# 1 + 0.75/3 for (A) and 1 + 0.75/1 for the others.
TOP_PAD_LINES = 0.75
TOP_VJUST_3LINE = 1 + TOP_PAD_LINES / 3
TOP_VJUST_1LINE = 1 + TOP_PAD_LINES
TOP_EXPAND = expansion(mult = c(0, 0.40))

# (A) Histogram of editing deltas ----
# Figure: histogram of the per-site editing-level shift introduced by WASP,
#         dEL = EL(WASP-passed) - EL(all reads), at each coverage floor.
#
# Reading it: POSITIVE dEL = editing was under-estimated without WASP, i.e. the
# direction reference bias predicts. The distribution is sharply peaked at 0
# with very small tails, meaning the WASP correction is negligible.
# The symmetric x range keeps the slight rightward skew visible.

dat = read_part(SD, "a_sites")
stats = read_part(SD, "a_stats")
stats[, lab := sprintf("Mean = %+.4f\nMedian = %+.4f\n%s sites",
                       avg, med, format(n, big.mark = ","))]

p.a = ggplot(dat, aes(dEL)) +
  geom_histogram(bins = 50, fill = unname(okabe_ito["blue"]), colour = "white",
                 linewidth = 0.1) +
  geom_text(data = stats, aes(x = -1, y = Inf, label = lab),
            hjust = 0, vjust = TOP_VJUST_3LINE, size = ANNOT_SIZE,
            lineheight = 1.1, inherit.aes = FALSE) +
  # Free y, deliberately: the evaluable set shrinks ~4x across the floors
  # (294,713 -> 167,647 -> 69,730 sites), so a shared axis would flatten the
  # cov >= 50 facet into a line and hide the SHAPE, which is the point.
  facet_wrap(~ cov_facet, nrow = 1, scales = "free_y") +
  # Linear y, deliberately: the point of this panel is that the correction is
  # negligible, which reads as a sharp spike at 0 with very small tails. A log
  # y-axis magnifies the tails and overstates the effect.
  scale_y_continuous(labels = shorten_axis, expand = TOP_EXPAND) +
  # coord_cartesian, not scale_x_continuous(limits=): scale limits drop rows
  # BEFORE binning, which silently discarded the some extreme sites.
  # Clipping at the coord level keeps every site in the histogram.
  coord_cartesian(xlim = c(-1, 1)) +
  labs(
    title = "Per-site editing level shifts",
    subtitle = "The correction WASP applies is negligible at every coverage floor",
    x = expression(Delta * " editing level (WASP-passed - all reads)"),
    y = "Number of sites"
  ) +
  theme_manuscript(base_size = BASE_SIZE, facets = TRUE) + grid_y

save_figure_sd(p.a, "WASP_hist_delta_editing_level", width = 9.6, height = 5.2, parts = c("a_sites", "a_stats"),
               grid = "y", script = script_id, src_stem = SD)

# (B) Read loss per allele ----
# Figure: fraction of reads removed by WASP at REDIportal A-to-I sites, split
#         by which allele the read carried (reference A vs edited G), at each
#         coverage floor.
#
# Reading it: this carries TWO claims at once. The bar HEIGHTS are the cost —
# WASP discards ~88% of the reads at the very sites being measured. The
# difference BETWEEN the bars is the reference-bias signature: reference-allele
# reads are removed at a marginally higher rate than edited-allele reads, which
# is the direction bias predicts, but the gap is sub-percentage-point and
# shrinks to nothing by cov >= 50. Removal here is close to indiscriminate
# rather than a surgical correction of an allele preference.
#
# This panel is expressed as RATES, not counts, on purpose: reference reads
# outnumber edited reads ~206:1, so on a count axis the edited bar would be
# half a percent of the reference bar's height and unreadable.
#
# The full 0-100% y axis is deliberate: a zoomed axis would dramatise a
# sub-percentage-point gap and hide that the dominant effect is bulk loss.

dat = read_part(SD, "b_alleles")
thr = read_part(SD, "b_gaps")

# Reference minus edited, per facet, with the CI from the two-proportion test.
gaps = merge(thr[, .(min_cov, gap_pp, gap_ci_lo_pp, gap_ci_hi_pp)],
             unique(dat[, .(min_cov, cov_facet)]), by = "min_cov")
# Gap only; its CI and p-value stay in report/<ds>/WASP_by_threshold.tsv rather
# than crowding the panel, which needs to carry one number per facet.
gaps[, lab := sprintf("Reference - Edited = %.2f pp", gap_pp)]

# The subtitle is DERIVED, never asserted. It used to end with a hardcoded
# "and it vanishes by cov >= 50", which is true of jzlxv3 and FALSE of jzlss2,
# where the gap is still 0.96 pp [0.89, 1.02] at the highest floor. A figure
# must not claim something its own panel contradicts, so whether the gap
# survives is read from the CI of the top threshold.
top_thr = thr[min_cov == max(min_cov)]
gap_vanishes = top_thr$gap_ci_lo_pp <= 0 & top_thr$gap_ci_hi_pp >= 0
sub.b = sprintf(
  "Bar height = the cost; the gap between bars = the bias, %s",
  if (isTRUE(gap_vanishes)) {
    sprintf("indistinguishable from zero by cov >= %d", top_thr$min_cov)
  } else {
    sprintf("still %.2f pp at cov >= %d", top_thr$gap_pp, top_thr$min_cov)
  })

p.b = ggplot(dat, aes(allele, removal_rate, fill = allele)) +
  geom_col(width = 0.6, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f%%\n(%s of\n%s reads)",
                                100 * removal_rate,
                                shorten_num(reads_removed),
                                shorten_num(reads_all))),
            vjust = -0.25, size = ANNOT_SIZE, lineheight = 1.1) +
  # At the top of the facet, not over the bars. The bars reach ~88%, so an
  # annotation placed mid-panel sits on top of them and needs a label box to
  # stay legible; the extra top expansion buys clear space instead.
  geom_text(data = gaps, aes(x = 1.5, y = Inf, label = lab),
            vjust = TOP_VJUST_1LINE, size = ANNOT_SIZE, inherit.aes = FALSE) +
  facet_wrap(~ cov_facet, nrow = 1) +
  scale_fill_manual_cvd(values = cmap_allele, guide = "none") +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = TOP_EXPAND) +
  labs(
    title = "Reads removed at REDIportal sites",
    subtitle = sub.b,
    x = NULL, y = "Percentage of allele's reads removed"
  ) +
  theme_manuscript(base_size = BASE_SIZE, facets = TRUE) + grid_y

save_figure_sd(p.b, "WASP_bar_removed_read_alleles", width = 9.6, height = 5.2, parts = c("b_alleles", "b_gaps"),
               grid = "y", script = script_id, src_stem = SD)

# (C) Site attrition ----
# Figure: REDIportal A-to-I sites adequately covered WITHOUT WASP, against the
#         sites still adequately covered after it, at each coverage floor.
#
# Reading it: this is the panel that decides the question. Sites carry no
# allele — every A-to-I site has both reference and edited reads — so this one
# is deliberately NOT split by allele; it counts sites, and (B) counts reads.
# WASP removes ~90% of the sites at which editing can be measured at all, to
# buy a correction of +0.0025 editing-level points. That is why WASP is
# reported as a diagnostic and kept out of the benchmark.

dat = read_part(SD, "c_sites")

# Shared y across facets, deliberately: the covered set itself shrinks as the
# floor rises (2.60M -> 1.77M -> 0.98M), and that shrinkage is part of the
# cost story. A free y would normalise it away.
p.c = ggplot(dat, aes(set, n_sites, fill = set)) +
  geom_col(width = 0.6, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = shorten_num(n_sites)),
            vjust = -0.35, size = ANNOT_SIZE) +
  # Same vjust and same top expansion as (B), so the gap between this text and
  # the bars matches that panel when the two are stacked.
  geom_text(data = unique(dat[, .(cov_facet, pct_sites_dropped)]),
            aes(x = 1.5, y = Inf, label = sprintf("-%.1f%% of sites", pct_sites_dropped)),
            vjust = TOP_VJUST_1LINE, size = ANNOT_SIZE, inherit.aes = FALSE) +
  facet_wrap(~ cov_facet, nrow = 1) +
  scale_fill_manual_cvd(values = cmap_attrition, guide = "none") +
  scale_y_continuous(labels = shorten_axis, expand = TOP_EXPAND) +
  labs(
    title = "REDIportal sites measurable before and after WASP",
    subtitle = "~90% of the sites where editing can be measured are lost to buy that correction",
    x = NULL, y = "Number of A-to-I sites"
  ) +
  theme_manuscript(base_size = BASE_SIZE, facets = TRUE) + grid_y

save_figure_sd(p.c, "WASP_bar_site_attrition", width = 9.6, height = 5.2, parts = "c_sites",
               grid = "y", script = script_id, src_stem = SD)

# (*) Combined figure panel ----
# Figure: (A), (B) and (C) stacked, the form the supplement figure takes.
#
# Stacked rather than side by side: all three share the same coverage-floor
# facets, so aligning them vertically lets a reader read one column straight
# down as "at this floor, the correction is X, the read cost is Y, the site
# cost is Z". ncol = 1 with align = "v" keeps those facet columns registered.

p.panel = cowplot::plot_grid(p.a, p.b, p.c, ncol = 1, align = "v", axis = "lr",
                             labels = "auto", label_size = 20)

save_figure_sd(p.panel, "WASP_panel_refbias", width = 9.6, height = 15.0,
               parts = c("a_sites", "a_stats", "b_alleles", "b_gaps", "c_sites"),
               grid = "y", script = script_id, src_stem = SD)
