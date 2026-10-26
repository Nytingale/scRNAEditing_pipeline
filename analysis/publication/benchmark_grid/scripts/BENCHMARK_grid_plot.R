#!/usr/bin/env Rscript
# ============================================================================
# BENCHMARK_grid_plot.R — every panel of the benchmark-grid deliverable (draws from source data only since
# 2026-10-25).
# ============================================================================
# Input : source_data/<DATASET>/BENCHMARK_grid.<part>.tsv (BENCHMARK_grid_source_data.R)
# Output: figures/<DATASET>/BENCHMARK_*.svg
# Run with:  bash scripts/run_BENCHMARK_figures.sh   or   Rscript BENCHMARK_grid_plot.R [DATASET]
#
# THE GRID. caller x bam_type x strand_mode x snp_subtract, over 20 samples.
# `sample` is the replicate dimension and is aggregated to a median with an IQR
# band; every other dimension is kept distinct by facet or aesthetic (checked by
# assert_no_dropped_combos() in BENCHMARK_grid_source_data.R).
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "BENCHMARK_plot_common.R"))
suppressPackageStartupMessages(library(cowplot))
script_id <- "publication/benchmark_grid/scripts/BENCHMARK_grid_plot.R"
SD <- "BENCHMARK_grid"
message("dataset: ", DATASET)

# ── (A) Final A-to-G fraction ───────────────────────────────────────────────
# Reading it: the headline specificity comparison. Real A-to-I editing appears
# as A->G, so a higher final A-to-G percentage means less of the surviving call
# set is other mismatch types. One box per (caller, bam_type); the box IS the
# across-sample distribution, which is why the source data keeps `sample` and
# holds the per-run values rather than a pre-computed median. Faceting on
# snp_subtract and strand_mode keeps those two settings separate: they were
# silently overplotted in the predecessor of this figure.
message("(A) BENCHMARK_box_ag_pct_final")
dat.a <- read_part(SD, "a_final")
p.a <- ggplot(dat.a, aes(x = bam_type, y = ag_pct_step7_pct, fill = caller)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, position = position_dodge(0.85)) +
  facet_grid(snp_subtract ~ strand_mode,
             labeller = labeller(snp_subtract = function(x) paste0("SNP subtract: ", x),
                                 strand_mode  = function(x) paste0("strand: ", x))) +
  scale_fill_manual_cvd(values = cmap_method, na.value = "grey50") +
  labs(x = "BAM preprocessing", y = "Final A-to-G sites (% of mismatches)",
       fill = "Caller") +
  theme_manuscript(facets = TRUE) + grid_y +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure_sd(p.a, "BENCHMARK_box_ag_pct_final", width = 11, height = 7.5, parts = "a_final", grid = "y",
               script = script_id, src_stem = SD)

# ── (B) Site attrition across filtering steps ───────────────────────────────
# Reading it: how many candidate sites each filter removes. The line is the
# across-sample median and the band the interquartile range, so a wide band
# means the filter behaves inconsistently between samples. y is log-scaled
# because step 1 and step 7 differ by orders of magnitude; the predecessor used
# free linear y per facet, which made the steep early drop invisible.
message("(B) BENCHMARK_line_site_attrition")
dat.b <- read_part(SD, "b_attrition")
attrition_panel <- function(d, ylab, logy) {
  p <- ggplot(d, aes(x = step_label, y = median, colour = caller, group = caller)) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = caller), alpha = 0.18,
                colour = NA, show.legend = FALSE) +
    geom_line(aes(linetype = snp_subtract), linewidth = 0.5) +
    geom_point(size = 1.1) +
    facet_grid(snp_subtract + strand_mode ~ bam_type) +
    scale_colour_manual_cvd(values = cmap_method, na.value = "grey50") +
    scale_fill_manual_cvd(values = cmap_method, na.value = "grey50") +
    labs(x = "Filtering step", y = ylab, colour = "Caller",
         linetype = "SNP subtraction") +
    theme_manuscript(facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (logy) p + scale_y_log10(labels = shorten_axis) else p + scale_y_continuous(labels = shorten_axis)
}
p.b <- attrition_panel(dat.b, "Surviving sites (median, IQR band)", logy = TRUE)
save_figure_sd(p.b, "BENCHMARK_line_site_attrition", width = 15, height = 10, parts = "b_attrition", grid = "y",
               script = script_id, src_stem = SD)

# ── (C) A-to-G fraction across filtering steps ──────────────────────────────
# Reading it: the counterpart to (B) — what each filter buys in specificity. (B)
# is the cost in sites, this is the benefit in A-to-G enrichment, so the pair
# shows whether a filter is worth its losses. Linear y, since the quantity is
# already a percentage.
message("(C) BENCHMARK_line_ag_pct_progression")
dat.c <- read_part(SD, "c_ag_pct")
p.c <- attrition_panel(dat.c, "A-to-G (%) (median, IQR band)", logy = FALSE)
save_figure_sd(p.c, "BENCHMARK_line_ag_pct_progression", width = 15, height = 10, parts = "c_ag_pct", grid = "y",
               script = script_id, src_stem = SD)

# ── (D) Strand mode: hybrid minus annotation ────────────────────────────────
# Reading it: the paired effect of hybrid strand assignment. Each point is one
# (caller, bam_type, snp_subtract) cell; the value is hybrid minus
# annotation-only A-to-G percentage, differenced WITHIN a run before any
# aggregation, so it is not contaminated by between-sample variation. The dashed
# zero line is "no difference": below it, hybrid makes specificity worse. Bars
# are the across-sample IQR of that paired difference.
# Only bcftools and gatk appear — the other callers emit no read-strand
# evidence, so hybrid would equal annotation and those runs are not produced.
message("(D) BENCHMARK_point_strand_mode_shift")
if (has_part(SD, "d_strand_shift")) {
  dat.d <- read_part(SD, "d_strand_shift")
  p.d <- ggplot(dat.d, aes(x = bam_type, y = d_ag_median, colour = caller)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(ymin = d_ag_q25, ymax = d_ag_q75), width = 0.25,
                  position = position_dodge(0.6), linewidth = 0.4) +
    geom_point(size = 2.2, position = position_dodge(0.6)) +
    facet_wrap(~snp_subtract, nrow = 1,
               labeller = labeller(snp_subtract = function(x) paste0("SNP subtract: ", x))) +
    scale_colour_manual_cvd(values = cmap_method, na.value = "grey50") +
    labs(x = "BAM preprocessing",
         y = "A-to-G % shift, hybrid - annotation\n(paired within run; median, IQR)",
         colour = "Caller") +
    theme_manuscript(facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_figure_sd(p.d, "BENCHMARK_point_strand_mode_shift", width = 8.5, height = 5.5, parts = "d_strand_shift",
                 grid = "y", script = script_id, src_stem = SD)
} else {
  message("    skipped: no d_strand_shift part (only one strand mode in this dataset)")
}

# ── (E) Adaptive filter windows ─────────────────────────────────────────────
# Reading it: which adaptive window each grid cell settled on — w1 for the step-5
# mismatch window, w2 for the step-6 cluster window — with the tile shaded by the
# step-6 A-to-G percentage those windows achieved. A tile cannot be dodged, so
# snp_subtract and strand_mode must be facets; two tiles in one cell is how the
# predecessor hid half of these values.
message("(E) BENCHMARK_heatmap_window")
if (has_part(SD, "e_windows")) {
  dat.e <- read_part(SD, "e_windows")
  p.e <- ggplot(dat.e, aes(x = caller, y = bam_type)) +
    geom_tile(aes(fill = ag_pct_step6_median), colour = "white", linewidth = 0.4) +
    geom_text(aes(label = paste0("w1=", w1_median, "\nw2=", w2_median)),
              size = 2.3, lineheight = 0.9) +
    facet_grid(snp_subtract ~ strand_mode,
               labeller = labeller(snp_subtract = function(x) paste0("SNP subtract: ", x),
                                   strand_mode  = function(x) paste0("strand: ", x))) +
    scale_fill_viridis_c(option = "D", name = "A-to-G %\n(step 6, median)") +
    labs(x = "Caller", y = "BAM preprocessing") +
    theme_manuscript(facets = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_figure_sd(p.e, "BENCHMARK_heatmap_window", width = 10, height = 8, parts = "e_windows", grid = "none",
                 script = script_id, src_stem = SD)
}

# ── (F) Where hybrid gets its strand call ───────────────────────────────────
# Reading it: in hybrid mode, how many final sites got their strand from read
# evidence versus the annotation fallback. It says how much of hybrid is actually
# read-driven, which is what (D)'s effect has to be attributed to. Stacked bars
# are why snp_subtract must be a facet here: with it dropped, the two runs
# stacked on top of each other and every bar was their SUM, not a value.
message("(F) BENCHMARK_stackedbar_strand_source")
if (has_part(SD, "f_strand_source")) {
  dat.f <- read_part(SD, "f_strand_source")
  p.f <- ggplot(dat.f, aes(x = caller, y = median, fill = source)) +
    geom_col(position = "stack", colour = "black", linewidth = 0.25) +
    facet_grid(snp_subtract ~ bam_type,
               labeller = labeller(snp_subtract = function(x) paste0("SNP subtract: ", x))) +
    scale_fill_manual_cvd(values = c(read = unname(cmap_strandMode[["hybrid"]]),
                                     annot = unname(cmap_strandMode[["annotation"]])),
                          labels = c(read = "Read evidence", annot = "Annotation"),
                          name = "Strand source") +
    scale_y_continuous(labels = shorten_axis) +
    labs(x = "Caller", y = "Final sites (median across samples)") +
    theme_manuscript(facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_figure_sd(p.f, "BENCHMARK_stackedbar_strand_source", width = 13, height = 6, parts = "f_strand_source",
                 grid = "y", script = script_id, src_stem = SD)
}

# ── (G) THE PUBLISHED FIGURE ────────────────────────────────────────────────
# Reading it: the two halves of one argument, so they must be read together —
# (a) what specificity each grid cell reaches, and (b) what the filters cost to
# get there. Judging either alone is how a cell that looks best on A-to-G but
# retains almost no sites gets mistaken for the winner. rel_heights match the
# heights the panels were saved at individually, so each keeps its own aspect
# ratio.
message("(G) BENCHMARK_panel_grid  [the published figure]")
p.panel <- plot_grid(p.a, p.c, ncol = 1, labels = "auto", rel_heights = c(7.5, 10))
save_figure_sd(p.panel, "BENCHMARK_panel_grid", width = 15, height = 17.5, parts = c("a_final", "c_ag_pct"),
               grid = "y", script = script_id, src_stem = SD)

message("\nAll panels written to ", FIG_DIR)
