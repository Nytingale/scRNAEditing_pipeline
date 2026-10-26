#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_gtdefinition_plot.R     (draws from source data only since 2026-10-25)
#
# Ground-truth threshold figures.
#   GTDEF_step1_quality          BQ/MQ 30 vs 20: substitution spectrum (no truth threshold)
#   GTDEF_step2_universe         dbSNP common: allele fractions + ADAR1-KO response (no truth threshold)
#   GTDEF_robustness             downstream thresholds: neighbour lattice + informative-KO bound
#   GTDEF_robust_r5_informative_ko  explainer for m
#   GTDEF_panel_sweeps           assembled: parameter sweeps barely move the ground truth or the benchmark
#   GTDEF_supp_gt_thresholds_v8  every ground-truth threshold on the v8 site set (+ GTDEF_v8_quality_universe)
#
# Input : source_data/pm_gt/GTDEF_gtdefinition.<part>.tsv (GTDEF_gtdefinition_source_data.R)
# Output: figures/pm_gt/GTDEF_step*.svg, GTDEF_robust*.svg, GTDEF_panel_sweeps.svg, GTDEF_supp_gt_thresholds_v8.svg,
#         GTDEF_v8_quality_universe.svg
#
# Colour meaning is fixed across every figure (plotting/plotting_functions.R):
#   cmap_gt_series  truth-set series (orange/red = depleted in KO, grey = not)
#   cmap_gt_overlap purple = overlap with the current choice's TRUE set
#   cmap_gt_redi    teal = % of TRUE sites in REDIportal
#
# Usage: Rscript GTDEF_gtdefinition_plot.R [STEP|all]    (STEP: step1_quality, step2_universe, robustness,
#        panel_sweeps, v8; default all)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

argv = commandArgs(trailingOnly = TRUE)
STEP_ARG = if (length(argv) >= 1) argv[1] else "all"
SCRIPT = "publication/gt_definition/scripts/GTDEF_gtdefinition_plot.R"
STEM = "GTDEF_gtdefinition"
# Parts are the report tables without the GTDEF_ prefix.
rp = function(stem) sd_file(SRC_DIR, STEM, sub("^GTDEF_", "", stem))
have = function(stem) file.exists(rp(stem))
rd = function(stem) as_tibble(read_part(STEM, sub("^GTDEF_", "", stem)))
PANEL_SWEEPS_NEEDS = c("step1_quality", "step2_universe", "robustness", "panel_sweeps")
want = function(st) STEP_ARG %in% c("all", st) || (STEP_ARG == "panel_sweeps" && st %in% PANEL_SWEEPS_NEEDS) ||
  (STEP_ARG == "v8" && st %in% c("step1_quality", "step2_universe"))
count_label = function(n) {
  out = formatC(n, format = "d", big.mark = ",")
  big = !is.na(n) & n >= 1000
  out[big] = vapply(n[big], shorten_num, character(1))
  out
}
BARW2 = 0.7


# (A) Shared panel builders  ----

series_lab = c(n_true = "Depleted in KO", n_true_powered = "Depleted in KO, powered",
               n_true_ci = "Depleted in KO, CI excludes threshold", n_fp = "Not depleted in KO")

mark_lab = function(lab, is_ref, status) if_else(is_ref, paste0(lab, "\n(", status, ")"), lab)

# One-dimensional sweep: truth counts / overlap with the current choice / % REDIportal
sweep_1d = function(d, xlab, title, value_lab, numeric_x = FALSE) {
  if (numeric_x) d = d %>% arrange(as.numeric(value))
  d = d %>% mutate(lab = mark_lab(value_lab(value), is_ref, status), lab = factor(lab, levels = unique(lab)))
  st = d$status[1]
  cnt = d %>% select(lab, n_true, n_true_powered, n_fp) %>%
    pivot_longer(-lab, names_to = "series", values_to = "n") %>%
    mutate(series = factor(series_lab[series], levels = series_lab[c("n_true", "n_true_powered", "n_fp")]))
  p1 = ggplot(cnt, aes(lab, n, fill = series)) +
    geom_col(position = position_dodge(width = BARW2 + 0.1), width = BARW2, colour = "grey20", linewidth = 0.2) +
    scale_fill_manual_cvd(values = cmap_gt_series[levels(cnt$series)], name = NULL) +
    scale_y_continuous(labels = shorten_axis) +
    labs(x = NULL, y = "Sites", title = title) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_manuscript() + grid_y + theme(legend.position = "top", axis.text.x = element_blank())
  ovl = c(jaccard_vs_ref = paste0("Jaccard with the ", st, " choice (TRUE sites)"),
          frac_ref_true_kept = paste0("Share of the ", st, " choice's TRUE sites retained"))
  ov = d %>% select(lab, all_of(names(ovl))) %>%
    pivot_longer(-lab, names_to = "metric", values_to = "v") %>%
    mutate(metric = factor(ovl[metric], levels = ovl))
  p2 = ggplot(ov, aes(lab, v, colour = metric, group = metric)) +
    geom_line(linewidth = 0.6) + geom_point(size = 2) +
    scale_colour_manual_cvd(values = setNames(unname(cmap_gt_overlap), ovl), name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "Overlap") + guides(colour = guide_legend(nrow = 1)) +
    theme_manuscript() + grid_y + theme(legend.position = "top", axis.text.x = element_blank())
  p3 = ggplot(d, aes(lab, pct_true_rediportal, group = 1)) +
    geom_line(colour = unname(cmap_gt_redi["rediportal"]), linewidth = 0.6) +
    geom_point(colour = unname(cmap_gt_redi["rediportal"]), size = 2) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = xlab, y = "TRUE in\nREDIportal (%)") +
    theme_manuscript() + grid_y
  cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr", rel_heights = c(1.6, 0.9, 0.8))
}

# Two-dimensional sweep heatmap; the current choice's tile is outlined.
heat_2d = function(d, x, y, fill, title, xlab, ylab, ramp = cmap_true_ramp) {
  ggplot(d, aes(.data[[x]], .data[[y]], fill = .data[[fill]])) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_tile(data = d %>% filter(is_ref), fill = NA, colour = "black", linewidth = 1.1) +
    geom_text(aes(label = count_label(.data[[fill]])), size = 3.2) +
    scale_fill_gradient(low = unname(ramp["low"]), high = unname(ramp["high"]), labels = shorten_axis, name = NULL) +
    labs(x = xlab, y = ylab, title = title,
         subtitle = paste0("Outlined: ", d$status[1], " choice")) +
    theme_manuscript(facets = TRUE) + theme(panel.grid = element_blank(), legend.position = "right")
}

save_step = function(p, stem, w, h, srcs, data = NULL, grid = "y")
  save_panel(p, file.path(FIG_DIR, stem), width = w, height = h, sources = srcs, grid = grid, script = SCRIPT,
             sidecar_dir = SRC_DIR)


# (Step 1) Pileup quality: substitution spectrum of every aligned base  ----
# No ground-truth threshold: every aligned base of merged NTC and ADAR1_KO.
if (want("step1_quality") && have("GTDEF_step1_spectrum")) {
  BS1 = 13
  q_lab = function(q) paste0("BQ/MQ ", q)
  Q_LEVELS = q_lab(c(20, 30))
  stopifnot(setequal(Q_LEVELS, names(cmap_bq_mq)))
  s_lab = c(merged_NTC_ctrl = "NTC", merged_ADAR1_KO = "ADAR1 KO")
  ty = rd("GTDEF_step1_spectrum")
  ty_levels = sort(unique(ty$type))
  ty = ty %>% mutate(q = factor(q_lab(bq_mq), levels = Q_LEVELS), type = factor(type, levels = ty_levels),
                     key = gsub(">", "", type), sample = factor(s_lab[sample], levels = s_lab))
  # (a) NTC: share of all mismatches by type, as in the GTDEF_bar_mismatch_spectrum plots (A>G highlighted)
  tn = ty %>% filter(sample == "NTC")
  p1a = ggplot(tn, aes(type, pct_of_mismatches / 100, fill = key)) +
    geom_col(width = 0.75, colour = "grey30", linewidth = 0.2) +
    geom_text(data = tn %>% filter(type == "A>G"), aes(label = sprintf("%.1f", pct_of_mismatches)),
              vjust = -0.4, size = 5.2) +
    scale_fill_manual(values = cmap_varType, guide = "none") +
    scale_y_continuous(labels = function(v) 100 * v, expand = expansion(mult = c(0, 0.15))) +
    facet_grid(. ~ q) +
    scale_x_discrete(labels = function(x) gsub(">", "", x)) +
    labs(x = "Mismatch type", y = "Percentage of NTC mismatches (%)") +
    theme_manuscript(base_size = BS1, facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # (b) ADAR dependence: ADAR1-KO rate / NTC rate (below 1 = lost without ADAR1)
  nk = rd("GTDEF_step1_ntc_vs_ko") %>%
    mutate(q = factor(q_lab(bq_mq), levels = Q_LEVELS), type = factor(type, levels = ty_levels))
  p1b = ggplot(nk, aes(type, ko_over_ntc, colour = q)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
    geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
    scale_colour_manual_cvd(values = cmap_bq_mq, name = NULL) +
    scale_y_continuous(limits = c(0.75, 1.05), breaks = seq(0.7, 1.1, 0.1)) +
    scale_x_discrete(labels = function(x) gsub(">", "", x)) +
    labs(x = "Mismatch type", y = "ADAR1 KO / NTC rate\n(per million bases)") +
    theme_manuscript(base_size = BS1) + grid_y + theme(legend.position = "top")
  sm = as_tibble(data.table::fread(rp("GTDEF_step1_summary"), integer64 = "double")) %>%
    mutate(sample = factor(s_lab[sample], levels = s_lab))
  kept = sm %>% select(sample, bq_mq, aligned_bases) %>% pivot_wider(names_from = bq_mq, values_from = aligned_bases) %>%
    mutate(pct = 100 * `30` / `20`)
  cap = paste0("Every aligned base, genome-wide (UMI-deduplicated merged BAMs).\nBQ/MQ 30 keeps ",
               paste0(sprintf("%.0f%%", kept$pct), " (", kept$sample, ")", collapse = ", "),
               " of the bases kept at BQ/MQ 20.")
  srcs = c(rp("GTDEF_step1_spectrum"), rp("GTDEF_step1_ntc_vs_ko"), rp("GTDEF_step1_summary"))
  fig = cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(cap, x = 0.02, hjust = 0, size = BS1 - 1, colour = "grey25"),
    p1a, p1b, ncol = 1, rel_heights = c(0.1, 0.9, 0.9), labels = c("", "a", "b"), label_size = 18)
  save_step(fig, "GTDEF_step1_quality", 11, 8.5, srcs, tn)
}


# (Step 2) Site universe: sites in dbSNP (dbSNP156 common) vs not  ----
# No ground-truth threshold: every candidate with NTC >= 10 reads.
if (want("step2_universe") && have("GTDEF_step2_af_hist")) {
  BS2 = 13
  sm = rd("GTDEF_step2_summary")
  grp_lab = sm %>% transmute(group, lab = paste0(group, " (", count_label(n_sites), ")"))
  relab_g = function(g) setNames(grp_lab$lab, grp_lab$group)[g]
  cols = setNames(unname(cmap_dbsnp[grp_lab$group]), grp_lab$lab)
  af = rd("GTDEF_step2_af_hist") %>% mutate(g = factor(relab_g(group), levels = grp_lab$lab))
  # (a) NTC allele fraction: germline SNPs sit near 50% / 100%, editing sites low
  p2a = ggplot(af, aes(100 * bin, 100 * frac, colour = g)) +
    geom_step(direction = "mid", linewidth = 0.9) +
    scale_colour_manual_cvd(values = cols, name = NULL) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0.01, 0.04))) +
    labs(x = "NTC A-to-G allele fraction (%)", y = "Share of sites (%)") +
    guides(colour = guide_legend(ncol = 1)) +
    theme_manuscript(base_size = BS2) + grid_y + theme(legend.position = "top")
  # (b) ADAR1-KO response: pooled KO / NTC A>G fraction per NTC allele-fraction bin
  # (1 = unchanged by ADAR1 loss); Wilson 95% CI of the pooled KO fraction
  bin_lab = function(b) { x = as.numeric(regmatches(b, gregexpr("[0-9.]+", b))[[1]]) * 100
    paste0(round(x[1]), "-", pmin(100, round(x[2])), "%") }
  kr = rd("GTDEF_step2_ko_response") %>%
    mutate(g = factor(relab_g(group), levels = grp_lab$lab), bl = map_chr(af_bin, bin_lab),
           bl = factor(bl, levels = unique(bl[order(ntc_edit)])))
  p2b = ggplot(kr, aes(bl, ko_over_ntc, colour = g, group = g)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
    geom_errorbar(aes(ymin = ko_lo / ntc_edit, ymax = ko_hi / ntc_edit), width = 0, linewidth = 0.6,
                  position = position_dodge(width = 0.3)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.3)) +
    annotate("text", x = 0.6, y = 1.05, label = "No change", hjust = 0, colour = "grey40", size = 5) +
    scale_colour_manual_cvd(values = cols, name = NULL, guide = "none") +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.25)) +
    labs(x = "NTC A-to-G allele fraction (bin)", y = "ADAR1 KO / NTC A-to-G fraction\n(pooled per bin)") +
    theme_manuscript(base_size = BS2) + grid_y
  cap = paste0("Every candidate site with NTC \u2265 ", sm$ntc_min_reads[1], " reads; no ground-truth threshold applied. ",
               "Sites in dbSNP = dbSNP156 common variants.\n",
               "A-to-G allele fraction = G reads / (A + G reads) at the site, on the transcript strand,
in the merged ",
               "UMI-deduplicated NTC (or ADAR1 KO) BAM, base and mapping quality \u2265 30.\n",
               "b: reads are pooled over all sites in each NTC allele-fraction bin, so no per-site KO coverage floor is needed.")
  srcs = c(rp("GTDEF_step2_summary"), rp("GTDEF_step2_af_hist"), rp("GTDEF_step2_ko_response"))
  fig = cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(cap, x = 0.02, hjust = 0, size = BS2 - 1, colour = "grey25"),
    cowplot::plot_grid(p2a, p2b, nrow = 1, rel_widths = c(1, 1.1), labels = c("a", "b"), label_size = 18,
                       align = "h", axis = "b"),
    ncol = 1, rel_heights = c(0.2, 1))
  save_step(fig, "GTDEF_step2_universe", 13, 7.2, srcs, af)
}


# (Robustness) How much does the ground-truth set change when one threshold moves?  ----
# Observed Jaccard values only; no stability cut-off is imposed. The informative
# KO threshold m is plotted on its own (GTDEF_robustness c), next to its bound (b).
if (want("robustness") && have("GTDEF_robust_centre_neighbours")) {
  # per-sample (per-teratoma) NTC coverage sits next to the reproducibility rule it feeds
  ax_lab = c(min_cov = "NTC coverage", min_edit_reads = "NTC edited reads", min_ctrl_edit = "NTC editing level",
             depletion = "Editing reduction in KO", min_cov_teratoma = "Per-sample NTC coverage",
             rule = "Reproducibility across samples", ko_min_expected = "Informative KO (m)")
  val_order = list(min_cov = as.character(8:12), min_cov_teratoma = as.character(8:12), ko_min_expected = as.character(1:5), min_edit_reads = as.character(2:4),
                   min_ctrl_edit = c("0.03", "0.04", "0.05", "0.06", "0.07"), depletion = c("0.3", "0.4", "0.5", "0.6", "0.7"),
                   rule = c("strict2", "evaluable", "strict3"))
  fmt_val = function(a, v) { a = rep_len(a, length(v)); v = as.character(v); case_when(a == "min_ctrl_edit" ~ paste0(round(100 * suppressWarnings(as.numeric(v))), "%"),
                                     a == "depletion" ~ paste0(round(100 * suppressWarnings(as.numeric(v))), "%"),
                                     a == "rule" & v == "strict2" ~ ">=2\nof 4", a == "rule" & v == "strict3" ~ ">=3\nof 4",
                                     a == "rule" ~ "all\nbut 1", TRUE ~ v) }
  order_values = function(df, vcol) df %>%
    mutate(v = as.character(.data[[vcol]]),
           v = if_else(axis %in% c("min_ctrl_edit", "depletion"), as.character(as.numeric(v)), v),
           ord = map2_int(axis, v, ~ match(.y, as.character(if (.x %in% c("min_ctrl_edit", "depletion"))
             as.numeric(val_order[[.x]]) else val_order[[.x]]))),
           lab = fmt_val(axis, v), axisf = factor(ax_lab[axis], levels = ax_lab)) %>%
    arrange(axisf, ord) %>% mutate(lab = factor(paste(axis, lab), levels = unique(paste(axis, lab))))
  relab = function(x) sub("^[a-z_]+ ", "", x)
  BS = 13                                   # base font size for the robustness figure
  GT_LAB = "Ground-truth sites"
  COL_SITES = unname(cmap_gt_series["Depleted in KO"]); COL_J = unname(cmap_gt_overlap["jaccard"])
  # One ggplot per row: facet per threshold, facet width proportional to its
  # number of values; x text horizontal (reproducibility labels on two lines).
  facet_row = function(p, ylab, x_text = TRUE)
    p + facet_grid(. ~ axisf, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 14)) +
      scale_x_discrete(labels = relab) + labs(x = NULL, y = ylab) +
      theme_manuscript(base_size = BS, facets = TRUE) + grid_y +
      theme(legend.position = "none", panel.spacing.x = unit(0.9, "lines")) +
      (if (!x_text) theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else theme())
  # ground-truth sites (top) and Jaccard with the centre (bottom) for a set of thresholds
  neighbour_rows = function(df) {
    ymax = max(df$n_neighbour) * 1.22
    top = facet_row(ggplot(df, aes(lab, n_neighbour, group = 1)) +
                      geom_line(colour = COL_SITES, linewidth = 0.8) +
                      geom_point(aes(size = centre), colour = COL_SITES) +
                      geom_point(data = df %>% filter(centre), colour = "black", shape = 21, size = 3.8, stroke = 1.1, fill = NA) +
                      geom_text(aes(label = count_label(n_neighbour)), vjust = -1.1, size = 3.6) +
                      scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 2.2), guide = "none") +
                      scale_y_continuous(labels = shorten_axis, limits = c(0, ymax), expand = expansion(mult = c(0, 0))),
                    GT_LAB, x_text = FALSE)
    bot = facet_row(ggplot(df, aes(lab, jaccard, group = 1)) +
                      geom_line(colour = COL_J, linewidth = 0.8) +
                      geom_point(colour = COL_J, size = 2.2) +
                      geom_text(data = df %>% filter(!centre), aes(label = sprintf("%.2f", jaccard)),
                                vjust = 1.9, size = 3.5, colour = "grey20") +
                      scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)),
                    "Jaccard with centre")
    cowplot::plot_grid(top, bot, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1.1))
  }
  title_row = function(txt)
    cowplot::ggdraw() + cowplot::draw_label(txt, x = 0.035, y = 0.5, hjust = 0, size = BS - 1, colour = "grey25")

  cn = rd("GTDEF_robust_centre_neighbours") %>% order_values("neighbour_value") %>% mutate(centre = step == 0)
  r1 = cowplot::plot_grid(
    title_row(paste0("Centre (circled): NTC >= 10 reads, >= 3 edited, >= 5% editing; >= 50% editing reduction in KO;\n",
                     "per-sample NTC >= 10 reads, edited in all but 1 evaluable sample (>= 2); m >= 3")),
    neighbour_rows(cn %>% filter(axis != "ko_min_expected")), ncol = 1, rel_heights = c(0.11, 1))
  rm_panel = neighbour_rows(cn %>% filter(axis == "ko_min_expected"))
  save_step(r1, "GTDEF_robust_r1_centre_neighbours", 15, 7, rp("GTDEF_robust_centre_neighbours"), cn)
}


# (Robustness, m) What the informative-KO threshold m means  ----
# m = KO coverage x NTC editing = edited KO reads expected with NO depletion.
# Reference lines (all panels): dashed = alpha 0.05 (m = 3.00), dotted = alpha 0.01 (m = 4.61).
if (want("robustness") && have("GTDEF_robust_m_theory")) {
  A_LT = c("alpha = 0.05" = "dashed", "alpha = 0.01" = "dotted")
  A_M = c("alpha = 0.05" = -log(0.05), "alpha = 0.01" = -log(0.01))
  ref = tibble(a = names(A_M), m = unname(A_M), p = exp(-unname(A_M)))
  th = rd("GTDEF_robust_m_theory") %>%
    mutate(lvl = factor(paste0(round(100 * ntc_edit), "%"), levels = names(cmap_ntc_edit_level)))
  # (a) chance of 0 edited KO reads with no depletion, exact binomial vs the exp(-m) bound
  m1 = ggplot(th, aes(m, p0_exact)) +
    geom_hline(data = ref, aes(yintercept = p, linetype = a), colour = "grey30") +
    geom_vline(data = ref, aes(xintercept = m, linetype = a), colour = "grey30") +
    geom_step(aes(colour = lvl), linewidth = 0.6) +
    geom_line(aes(y = p0_bound), colour = "black", linewidth = 0.9) +
    annotate("text", x = 0.9, y = exp(-0.9) * 1.6, label = "exp(-m)", hjust = 0, size = 3) +
    scale_colour_manual_cvd(values = cmap_ntc_edit_level, name = "NTC editing p") +
    scale_linetype_manual(values = A_LT, name = NULL) +
    scale_y_log10(breaks = c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002)) +
    scale_x_continuous(breaks = 0:6) +
    labs(x = "m = KO coverage x NTC editing", y = "P(0 edited KO reads | no depletion)",
         title = "a - m bounds the chance that a KO arm shows 0 edited reads by luck",
         subtitle = "Steps: exact (1 - p)^n at the smallest KO coverage n with n x p >= m
Black: exp(-m), an upper bound for every p") +
    theme_manuscript() + grid_both + theme(legend.position = "right")
  # (b) KO coverage a site needs, by its NTC editing level
  kc = rd("GTDEF_robust_ko_cov_needed") %>%
    mutate(a = if_else(abs(alpha - 0.05) < 1e-9, "alpha = 0.05", "alpha = 0.01")) %>%
    pivot_longer(c(ko_cov_m_rule, ko_cov_exact), names_to = "rule", values_to = "ko_cov") %>%
    mutate(rule = if_else(rule == "ko_cov_m_rule", "m rule: n >= m / p", "exact: (1 - p)^n <= alpha"))
  m2 = ggplot(kc, aes(100 * ntc_edit, ko_cov, linetype = a, colour = rule)) +
    geom_step(linewidth = 0.6) +
    scale_linetype_manual(values = A_LT, name = NULL) +
    scale_colour_manual_cvd(values = c("m rule: n >= m / p" = "black", "exact: (1 - p)^n <= alpha" = "grey60"), name = NULL) +
    scale_x_log10(breaks = c(3, 5, 10, 20, 50, 80)) + scale_y_log10(breaks = c(3, 5, 10, 20, 50, 100, 200)) +
    labs(x = "NTC editing p (%, log)", y = "KO coverage needed (log)",
         title = "b - KO coverage a site needs to be informative",
         subtitle = "e.g. p = 10%: 30 KO reads (alpha 0.05) or 47 (alpha 0.01)") +
    theme_manuscript() + grid_both + theme(legend.position = "right")
  srcs = c(rp("GTDEF_robust_m_theory"), rp("GTDEF_robust_ko_cov_needed"))
  panels = list(m1, m2)
  if (have("GTDEF_robust_m_curve")) {
    mc = rd("GTDEF_robust_m_curve")
    mcl = mc %>% pivot_longer(c(n_gated, n_true, n_true_powered), names_to = "s", values_to = "n") %>%
      mutate(s = factor(c(n_gated = "Gated, informative KO", n_true = "Depleted in KO",
                          n_true_powered = "Depleted in KO, powered")[s],
                        levels = c("Gated, informative KO", "Depleted in KO", "Depleted in KO, powered")))
    at = mc %>% filter(abs(ko_min_expected - A_M[1]) < 1e-3 | abs(ko_min_expected - A_M[2]) < 1e-3) %>%
      mutate(lab = paste0("m = ", sprintf("%.2f", ko_min_expected), ": ", count_label(n_true), " sites"))
    m3 = ggplot(mcl, aes(ko_min_expected, n, colour = s)) +
      geom_vline(data = ref, aes(xintercept = m, linetype = a), colour = "grey30") +
      geom_line(linewidth = 0.8) +
      geom_point(data = at, aes(ko_min_expected, n_true), inherit.aes = FALSE,
                 colour = unname(cmap_gt_series["Depleted in KO"]), size = 2.4) +
      geom_text(data = at, aes(ko_min_expected, n_true, label = lab), inherit.aes = FALSE,
                hjust = -0.08, vjust = -0.8, size = 2.8) +
      scale_colour_manual_cvd(values = c("Gated, informative KO" = unname(cmap_grey_ramp["high"]),
                                         "Depleted in KO" = unname(cmap_gt_series["Depleted in KO"]),
                                         "Depleted in KO, powered" = unname(cmap_gt_series["Depleted in KO, powered"])),
                              name = NULL) +
      scale_linetype_manual(values = A_LT, name = NULL) + guides(linetype = guide_legend(order = 1)) +
      scale_x_continuous(breaks = 0:6) + scale_y_log10(labels = shorten_axis) +
      labs(x = "m = KO coverage x NTC editing", y = "Sites (log)",
           title = "c - sites kept as m increases (all other thresholds at the centre)",
           subtitle = "m = 0: any KO coverage; NTC >= 10 reads, >= 3 edited, >= 5%, evaluable, d = 50%") +
      theme_manuscript() + grid_both + theme(legend.position = "right")
    srcs = c(srcs, rp("GTDEF_robust_m_curve")); panels = c(panels, list(m3))
  }
  save_step(cowplot::plot_grid(plotlist = panels, ncol = 1, align = "v", axis = "lr"),
            "GTDEF_robust_r5_informative_ko", 10.5, 4.4 * length(panels), srcs)
}


# (Robustness) Combined figure: a = thresholds; b = the m bound; c = the m threshold  ----
if (want("robustness") && exists("r1") && exists("rm_panel")) {
  # (b) a single curve: P(0 edited KO reads) under the no-depletion null, Poisson exp(-m)
  mb = tibble(m = seq(0, 6, 0.02)) %>% mutate(p0 = exp(-m))
  mref = tibble(m = -log(c(0.05, 0.01)), p0 = c(0.05, 0.01),
                lab = c("m = 3: P = 0.05 (chosen)", "m = 4.6: P = 0.01"))
  Y0 = exp(-6)
  mb_panel = ggplot(mb, aes(m, p0)) +
    geom_segment(data = mref, aes(x = m, xend = m, y = Y0, yend = p0), linetype = 2, colour = "grey45") +
    geom_segment(data = mref, aes(x = 0, xend = m, y = p0, yend = p0), linetype = 2, colour = "grey45") +
    geom_line(linewidth = 1.1, colour = "black") +
    geom_point(data = mref, size = 3, colour = "black") +
    geom_text(data = mref, aes(label = lab), hjust = -0.08, vjust = -0.7, size = 4.2) +
    scale_y_log10(breaks = c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005), limits = c(Y0, 1), expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(breaks = 0:6, limits = c(0, 6.6), expand = expansion(mult = c(0, 0))) +
    labs(x = "m = expected edited KO reads if editing is unchanged\n(KO coverage x NTC editing level)",
         y = "P(0 edited KO reads)",
         title = NULL) +
    theme_manuscript(base_size = BS) + grid_both +
    theme(plot.margin = margin(5, 10, 5, 5))
  mb_panel = cowplot::plot_grid(title_row("Chance of seeing no edited KO reads when editing is unchanged"), mb_panel,
                                ncol = 1, rel_heights = c(0.09, 1))
  mc_panel = cowplot::plot_grid(title_row("Ground-truth set as the informative-KO threshold m moves"), rm_panel,
                                ncol = 1, rel_heights = c(0.09, 1))
  save_step(cowplot::plot_grid(r1, cowplot::plot_grid(mb_panel, mc_panel, nrow = 1, rel_widths = c(1.15, 1),
                                                      labels = c("b", "c"), label_size = 18),
                               ncol = 1, rel_heights = c(1, 0.95), labels = c("a", ""), label_size = 18),
            "GTDEF_robustness", 15, 12.5,
            rp("GTDEF_robust_centre_neighbours"), mb)
}


# (Panel) Parameter sweeps barely move the ground truth or the benchmark  ----
#   a-b  BQ/MQ 30 vs 20 (step 1)          c-d  dbSNP common sites (step 2)
#   e    one threshold at a time: ground-truth sites + Jaccard with the centre (m excluded)
#   f    the informative-KO bound: P(0 edited KO reads) = exp(-m) (GTDEF_robustness b)
#   g    best configuration per caller, F1 vs depletion threshold, per tier
#   h    Spearman of configuration F1 between depletion thresholds, per tier
# The pipeline's window filters are benchmarking filters, not ground-truth criteria,
# so truth sets with those filters applied are not shown here.
# g-h: parts sweep_caller_best / sweep_rank_stability (gt_sweep_infneg_v6.R).
if (want("panel_sweeps") && all(c("p1a", "p1b", "p2a", "p2b", "cn", "mb_panel") %in% ls())) {
  BSP = 13
  tier_lab = c(arm = "Arm", sample = "Sample", celltype = "Cell type")
  t_lab = function(t) paste0(sub("^t", "", t), "%")

  cb = rd("GTDEF_sweep_caller_best") %>% distinct(tier, t, caller, F_real) %>%
    mutate(tier = factor(tier_lab[tier], levels = tier_lab), d = as.numeric(sub("^t", "", t)))
  ph = ggplot(cb, aes(d, F_real, colour = caller, group = caller)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    facet_wrap(~ tier, nrow = 1) +
    scale_colour_manual(values = cmap_method[sort(unique(cb$caller))], labels = label_caller, name = NULL) +
    scale_x_continuous(breaks = c(25, 50, 75), labels = function(x) paste0(x, "%"), expand = expansion(mult = 0.12)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Editing reduction in KO (positives)", y = "F1, best configuration\nper caller") +
    theme_manuscript(base_size = BSP, facets = TRUE) + grid_y + theme(legend.position = "bottom")
  rs = rd("GTDEF_sweep_rank_stability") %>%
    mutate(tier = factor(tier_lab[tier], levels = tier_lab), pair = paste(t_lab(t_a), "vs", t_lab(t_b)))
  p_i = ggplot(rs, aes(pair, rho_f1)) +
    geom_col(fill = unname(cmap_gt_overlap["jaccard"]), width = 0.65) +
    geom_text(aes(label = sprintf("%.3f", rho_f1)), vjust = -0.5, size = 3.6) +
    facet_wrap(~ tier, nrow = 1) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
    labs(x = "Depletion thresholds compared", y = "Spearman \u03C1 of configuration F1") +
    theme_manuscript(base_size = BSP, facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  row_lab = function(txt) cowplot::ggdraw() +
    cowplot::draw_label(txt, x = 0.01, hjust = 0, size = BSP + 1, fontface = "bold", colour = "grey15")
  LS = 18

  r_a = cowplot::plot_grid(p1a, p1b, nrow = 1, labels = c("a", "b"), label_size = LS, rel_widths = c(1.25, 1),
                           align = "h", axis = "b")
  r_c = cowplot::plot_grid(p2a, p2b, nrow = 1, labels = c("c", "d"), label_size = LS, rel_widths = c(1, 1.1),
                           align = "h", axis = "b")
  r_e = cowplot::plot_grid(neighbour_rows(cn %>% filter(axis != "ko_min_expected")), mb_panel, nrow = 1,
                           labels = c("e", "f"), label_size = LS, rel_widths = c(2.5, 1.2))
  r_h = cowplot::plot_grid(ph, p_i, nrow = 1, labels = c("g", "h"), label_size = LS, rel_widths = c(1.2, 1),
                           align = "h", axis = "b")
  fig = cowplot::plot_grid(
    row_lab("Pileup quality: base and mapping quality 30 vs 20 (every aligned base)"), r_a,
    row_lab("Site universe: dbSNP common sites behave as germline variants, not editing"), r_c,
    row_lab(paste0("One threshold at a time around the centre: NTC >= 10 reads, >= 3 edited, >= 5% editing, ",
                   ">= 50% reduction in KO, per-sample NTC >= 10, all but 1 sample, m >= 3; and the m bound")), r_e,
    row_lab("Benchmark: method ranking across depletion thresholds (informative negatives)"), r_h,
    ncol = 1, rel_heights = c(0.05, 0.9, 0.05, 0.85, 0.05, 1, 0.05, 0.8))
  srcs = c(rp("GTDEF_step1_spectrum"), rp("GTDEF_step1_ntc_vs_ko"), rp("GTDEF_step1_summary"),
           rp("GTDEF_step2_summary"), rp("GTDEF_step2_af_hist"), rp("GTDEF_step2_ko_response"),
           rp("GTDEF_robust_centre_neighbours"),
           rp("GTDEF_sweep_caller_best"), rp("GTDEF_sweep_rank_stability"))
  save_step(fig, "GTDEF_panel_sweeps", 20, 24.5, srcs, cn)
}


# (v8) Supplementary figure: every ground-truth threshold, on the v8 site set  ----
#   a-b  BQ/MQ 30 vs 20 (every aligned base)      c-d  dbSNP common sites (every candidate)
#   e    one threshold at a time around the centre, both gates (gt_robustness_v8.R)
#   f    P(0 edited KO reads) = exp(-m)           g    KO coverage a site needs to be informative
#   h    expected gate: sites vs m                i    depth gate: sites vs KO coverage
#   j    best configuration per caller, F1 vs depletion threshold
#   k    Spearman of configuration F1 between depletion thresholds
# Replaces GTDEF_panel_sweeps / GTDEF_robustness / GTDEF_robust_r1_centre_neighbours /
# GTDEF_robust_r5_informative_ko (v6 site set) and carries GTDEF_step1_quality / GTDEF_step2_universe.
if (STEP_ARG %in% c("all", "v8") && have("GTDEF_robust_centre_neighbours_v8") &&
    all(c("p1a", "p1b", "p2a", "p2b") %in% ls())) {
  BSP = 13; LS = 18
  GATE_COL = c(expected = "#0072B2", depth = "#009E73")
  GATE_LAB = V8_GATE_LAB
  gate_scale = function(aes = "colour")
    scale_discrete_manual(aes, values = setNames(GATE_COL, GATE_LAB), name = NULL)
  ax8 = c(min_cov = "NTC coverage", min_edit_reads = "NTC edited reads", min_ctrl_edit = "NTC editing level",
          min_cov_teratoma = "Per-sample NTC coverage", rule = "Reproducibility across samples",
          depletion = "Editing reduction in KO", ko_min_expected = "Expected edited KO reads (m)",
          ko_min_cov = "KO coverage")
  fmt8 = function(a, v) case_when(a %in% c("min_ctrl_edit", "depletion") ~ paste0(round(100 * suppressWarnings(as.numeric(v))), "%"),
                                  a == "rule" & v == "strict2" ~ "\u22652\nof 4", a == "rule" & v == "strict3" ~ "\u22653\nof 4",
                                  a == "rule" ~ "all\nbut 1", TRUE ~ as.character(v))
  cn8 = rd("GTDEF_robust_centre_neighbours_v8") %>%
    mutate(lab = fmt8(axis, neighbour_value), axisf = factor(ax8[axis], levels = ax8), centre = step == 0,
           g = factor(GATE_LAB[gate], levels = GATE_LAB)) %>%
    arrange(axisf, step) %>% mutate(lab = factor(paste(axis, step, lab), levels = unique(paste(axis, step, lab))))
  relab8 = function(x) sub("^[a-z_]+ -?[0-9]+ ", "", x)
  frow = function(p, ylab, x_text = TRUE)
    p + facet_grid(. ~ axisf, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 14)) +
      scale_x_discrete(labels = relab8) + labs(x = NULL, y = ylab) + gate_scale() +
      theme_manuscript(base_size = BSP, facets = TRUE) + grid_y +
      theme(legend.position = "none", panel.spacing.x = unit(0.9, "lines")) +
      (if (!x_text) theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else theme())
  e_top = frow(ggplot(cn8, aes(lab, n_neighbour, colour = g, group = g)) +
                 geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
                 geom_point(data = cn8 %>% filter(centre), shape = 21, size = 3.8, stroke = 1.1, fill = NA, colour = "black") +
                 geom_text(aes(label = count_label(n_neighbour), vjust = if_else(gate == "depth", -1, 1.9)),
                           size = 3.6, show.legend = FALSE) +
                 scale_y_continuous(labels = shorten_axis, limits = c(0, max(cn8$n_neighbour) * 1.15)),
               "Ground-truth sites\n(\u2265 50% editing\nreduction in KO)", x_text = FALSE) +
    theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1))
  e_bot = frow(ggplot(cn8, aes(lab, jaccard, colour = g, group = g)) +
                 geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
                 geom_text(data = cn8 %>% filter(!centre), aes(label = sprintf("%.2f", jaccard),
                           vjust = if_else(gate == "depth", -0.9, 1.9)), size = 3.4, show.legend = FALSE) +
                 geom_hline(yintercept = 1, colour = "grey80", linewidth = 0.3) +
                 scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)),
               "Jaccard with centre")
  p_e = cowplot::plot_grid(e_top, e_bot, ncol = 1, align = "v", axis = "lr", rel_heights = c(1.15, 1))

  # f: the no-depletion null, Poisson exp(-m)
  mb = tibble(m = seq(0, 6, 0.02)) %>% mutate(p0 = exp(-m))
  mref = tibble(m = -log(c(0.05, 0.01)), p0 = c(0.05, 0.01), lab = c("m = 3: P = 0.05\n(chosen)", "m = 4.6: P = 0.01"),
                hj = c(-0.08, 1.04), vj = c(-0.25, 1.5))
  Y0 = exp(-6)
  p_f = ggplot(mb, aes(m, p0)) +
    geom_segment(data = mref, aes(x = m, xend = m, y = Y0, yend = p0), linetype = 2, colour = "grey45") +
    geom_segment(data = mref, aes(x = 0, xend = m, y = p0, yend = p0), linetype = 2, colour = "grey45") +
    geom_line(linewidth = 1.1) + geom_point(data = mref, size = 3) +
    geom_text(data = mref, aes(label = lab, hjust = hj, vjust = vj), size = 5, lineheight = 0.9) +
    scale_y_log10(breaks = c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005), limits = c(Y0, 1)) +
    scale_x_continuous(breaks = 0:6, limits = c(0, 7)) +
    labs(x = "m = expected edited KO reads\n(KO coverage \u00D7 NTC editing level)", y = "P(0 edited KO reads | editing unchanged)") +
    theme_manuscript(base_size = BSP) + grid_both
  # g: KO coverage needed (theory; GTDEF_robust_r5_informative_ko b)
  A_LT = c("\u03B1 = 0.05" = "dashed", "\u03B1 = 0.01" = "dotted")
  kc = rd("GTDEF_robust_ko_cov_needed") %>%
    mutate(a = if_else(abs(alpha - 0.05) < 1e-9, "\u03B1 = 0.05", "\u03B1 = 0.01")) %>%
    pivot_longer(c(ko_cov_m_rule, ko_cov_exact), names_to = "rule", values_to = "ko_cov") %>%
    mutate(rule = if_else(rule == "ko_cov_m_rule", "m rule: n \u2265 m / p", "exact: (1 \u2212 p)\u207F \u2264 \u03B1"))
  p_g = ggplot(kc, aes(100 * ntc_edit, ko_cov, linetype = a, colour = rule)) +
    geom_hline(yintercept = 10, colour = GATE_COL[["depth"]], linewidth = 0.5) +
    geom_step(linewidth = 0.6) +
    scale_linetype_manual(values = A_LT, name = NULL) +
    scale_colour_manual(values = c("m rule: n \u2265 m / p" = "black", "exact: (1 \u2212 p)\u207F \u2264 \u03B1" = "grey60"), name = NULL) +
    scale_x_log10(breaks = c(3, 5, 10, 20, 50, 80)) + scale_y_log10(breaks = c(3, 5, 10, 20, 50, 100, 200)) +
    labs(x = "NTC editing level p (%, log)", y = "KO coverage needed (log)") +
    guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
    theme_manuscript(base_size = BSP) + grid_both + theme(legend.position = "bottom", legend.box = "vertical")
  # h-i: sites as the informative-KO threshold of each gate moves (v8 site set)
  cs_lab = c(n_gated = "All sites passing the coverage gate", n_true_t25 = "\u2265 25% editing reduction in KO",
             n_true = "\u2265 50% editing reduction in KO", n_true_t75 = "\u2265 75% editing reduction in KO")
  cs_col = c(unname(cmap_grey_ramp["high"]), "#FDAE6B", "#E6550D", "#A63603")
  cs_lt = c("solid", "solid", "solid", "solid")
  curve_panel = function(d, x, chosen, xlab, ttl, x_step, cols = cs_col, sl = cs_lab) {
    keep = intersect(names(sl), names(d))
    dl = d %>% pivot_longer(all_of(keep), names_to = "s", values_to = "n") %>%
      mutate(s = factor(sl[s], levels = sl[keep]))
    at = d %>% filter(abs(.data[[x]] - chosen) < 1e-6)
    ggplot(dl, aes(.data[[x]], log10(n), colour = s, linetype = s)) +
      geom_vline(xintercept = chosen, linetype = 2, colour = "grey45") +
      geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
      scale_linetype_manual(values = setNames(cs_lt, sl)[sl[keep]], name = NULL) +
      scale_colour_manual(values = setNames(cols, sl), name = NULL) +
      scale_x_continuous(breaks = seq(0, max(d[[x]]), x_step)) +
      labs(x = xlab, y = expression(Log[10]~"Number of sites"), title = ttl) +
      guides(colour = guide_legend(ncol = 1)) +
      theme_manuscript(base_size = BSP) + grid_both + theme(legend.position = "bottom", plot.title = element_text(size = BSP))
  }
  cs_lab0 = c(n_true_t0 = "\u2265 0% editing reduction in KO", cs_lab[-1])
  CS_EXP = c("black", "#99C7E0", "#0072B2", "#003F62"); CS_DEP = c("black", "#99D8C7", "#009E73", "#00573F")
  leg_bl = theme(legend.position = "inside", legend.position.inside = c(0.01, 0.01), legend.justification = c(0, 0),
                 legend.background = element_blank(), legend.key.size = unit(0.9, "lines"),
                 legend.text = element_text(size = BSP - 1))
  leg_room = function(d) expand_limits(y = log10(min(d$n_true_t75)) - 0.55)
  mcv = rd("GTDEF_robust_m_curve_v8"); kcv = rd("GTDEF_robust_ko_cov_curve_v8")
  p_h = curve_panel(mcv, "ko_min_expected", 3, "m = expected edited KO reads", V8_GATE_LAB2[["expected"]], 1,
                    cols = CS_EXP, sl = cs_lab0) + leg_bl + leg_room(mcv)
  p_i = curve_panel(kcv, "ko_min_cov", 10, "KO coverage (reads)", V8_GATE_LAB2[["depth"]], 5,
                    cols = CS_DEP, sl = cs_lab0) + leg_bl + leg_room(kcv)

  # j-k: benchmark sensitivity to the depletion threshold (v8 scoring, call state mismatch)
  tier_lab = c(arm = "Arm", sample = "Sample", celltype = "Cell type")
  t_lab = function(t) paste0(sub("^t", "", t), "%")
  cb = rd("GTDEF_sweep_caller_best_v8") %>%
    mutate(tier = factor(tier_lab[tier], levels = tier_lab), d = as.numeric(sub("^t", "", t)),
           g = factor(GATE_LAB[gate], levels = GATE_LAB))
  p_j = ggplot(cb, aes(d, f1, colour = caller, group = caller)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    facet_grid(g ~ tier, labeller = labeller(g = label_wrap_gen(14))) +
    scale_colour_manual(values = cmap_method[sort(unique(cb$caller))], labels = label_caller, name = NULL) +
    scale_x_continuous(breaks = c(25, 50, 75), labels = function(x) paste0(x, "%"), expand = expansion(mult = 0.12)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Editing reduction in KO", y = "F1, best configuration per caller") +
    theme_manuscript(base_size = BSP, facets = TRUE) + grid_y + theme(legend.position = "bottom")
  rs = rd("GTDEF_sweep_rank_stability_v8") %>%
    mutate(tier = factor(tier_lab[tier], levels = tier_lab), pair = paste(t_lab(t_a), "vs", t_lab(t_b)),
           g = factor(GATE_LAB[gate], levels = GATE_LAB))
  p_k = ggplot(rs, aes(pair, rho_f1, fill = g)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", rho_f1)), position = position_dodge(width = 0.75), vjust = -0.4, size = 3) +
    facet_wrap(~ tier, nrow = 1) + gate_scale("fill") +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
    labs(x = "Editing-reduction thresholds compared", y = "Spearman's \u03C1 of configuration F1") +
    theme_manuscript(base_size = BSP, facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

  row_lab = function(txt) cowplot::ggdraw() +
    cowplot::draw_label(txt, x = 0.01, hjust = 0, size = BSP + 1, fontface = "bold", colour = "grey15")
  cl8 = rd("GTDEF_robust_cleaning_v8")
  r_a = cowplot::plot_grid(p1a, p1b, nrow = 1, labels = c("a", "b"), label_size = LS, rel_widths = c(1.25, 1), align = "h", axis = "b")
  r_c = cowplot::plot_grid(p2a, p2b, nrow = 1, labels = c("c", "d"), label_size = LS, rel_widths = c(1, 1.1), align = "h", axis = "b")
  r_e = cowplot::plot_grid(p_e, labels = "e", label_size = LS)
  r_f = cowplot::plot_grid(p_f, p_g, p_h, p_i, nrow = 1, labels = c("f", "g", "h", "i"), label_size = LS,
                           align = "h", axis = "t")
  r_j = cowplot::plot_grid(p_j, p_k, nrow = 1, labels = c("j", "k"), label_size = LS, rel_widths = c(1.15, 1),
                           align = "h", axis = "b")
  fig = cowplot::plot_grid(
    row_lab("Pileup quality: base and mapping quality 30 vs 20 (every aligned base; no truth threshold)"), r_a,
    row_lab("Site universe: dbSNP common sites behave as germline variants, not editing (every candidate; no truth threshold)"), r_c,
    row_lab(paste0("One threshold at a time around the centre (RNA A-to-G; multiallelic sites, PGP1 SNVs and sites within ",
                   cl8$w1[1], " bp of other mismatch types removed)\nCentre: NTC \u2265 10 reads, \u2265 3 edited, ",
                   "\u2265 5% editing; per-sample NTC \u2265 10, edited in all but 1 evaluable sample; \u2265 50% reduction ",
                   "in KO; power-gated (m \u2265 3 expected edited KO reads) or depth-gated (\u2265 10 KO reads). Circled = centre")), r_e,
    row_lab("Informative KO arm: the no-depletion null and the sites each gate keeps"), r_f,
    row_lab("Benchmark: method ranking across editing-reduction thresholds (call sets: sites near other mismatch types removed)"), r_j,
    ncol = 1, rel_heights = c(0.04, 0.8, 0.04, 0.75, 0.07, 1.5, 0.04, 0.85, 0.04, 1.05))
  srcs = c(rp("GTDEF_step1_spectrum"), rp("GTDEF_step1_ntc_vs_ko"), rp("GTDEF_step1_summary"),
           rp("GTDEF_step2_summary"), rp("GTDEF_step2_af_hist"), rp("GTDEF_step2_ko_response"),
           rp("GTDEF_robust_centre_neighbours_v8"), rp("GTDEF_robust_cleaning_v8"),
           rp("GTDEF_robust_ko_cov_needed"), rp("GTDEF_robust_m_curve_v8"), rp("GTDEF_robust_ko_cov_curve_v8"),
           rp("GTDEF_sweep_caller_best_v8"), rp("GTDEF_sweep_rank_stability_v8"))
  save_step(fig, "GTDEF_supp_gt_thresholds_v8", 22, 32, srcs, cn8)
  # panels a-d alone (review response: pileup quality and site universe)
  fig_ad = cowplot::plot_grid(r_a, r_c, ncol = 1, rel_heights = c(1, 0.95))
  save_step(fig_ad, "GTDEF_v8_quality_universe", 12, 10,
            c(rp("GTDEF_step1_spectrum"), rp("GTDEF_step1_ntc_vs_ko"), rp("GTDEF_step1_summary"),
              rp("GTDEF_step2_summary"), rp("GTDEF_step2_af_hist"), rp("GTDEF_step2_ko_response")))
  # Review response figures 1 and 2 (GTDEF_v8_informative_ko, GTDEF_v8_universe_celltype_robustness): since
  # 2026-10-25 in GTDEF_v8_informative_ko_{source_data,plot}.R and GTDEF_v8_universe_celltype_robustness_*.R.
}


message("\nfigures (", STEP_ARG, ") -> ", FIG_DIR)
