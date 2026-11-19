#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_truth_spectrum12_plot.R   (user 2026-11-12: every ground truth gets a mismatch spectrum)
#
# Draws the source data of GTDEF_truth_spectrum12_source_data.R:
#   GTDEF_bar_truth_spectrum12_<dataset>_<truth>   one figure per ground truth: % of KO-depleted gate sites per
#                                                  mismatch type; rows = gate x before / after cleaning,
#                                                  columns = editing reduction in KO >= 25 / 50 / 75%
#   GTDEF_bar_truth_spectrum12_<dataset>_truths    the five truths of a dataset after cleaning, t = 50%
# The v9 positives are the A>G and T>C rows (read orientation) after cleaning.
#
# Input : source_data/pm_gt/GTDEF_truth_spectrum12.spectrum.tsv
# Output: figures/pm_gt/GTDEF_bar_truth_spectrum12_*.{svg,pdf}
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
SCRIPT <- "publication/gt_definition/scripts/GTDEF_truth_spectrum12_plot.R"
STEM <- "GTDEF_truth_spectrum12"
BS <- 12
TYPES <- all_vartypes()
GATES <- c("Expected-count gate", "Depth gate")
STAGES <- c("Before filtering", "After filtering")

sp <- read_part(STEM, "spectrum")
sp[, `:=`(type = factor(type, levels = TYPES),
          t_lab = factor(sprintf("Editing reduction in KO \u2265 %s%%", sub("^t", "", t)),
                         levels = sprintf("Editing reduction in KO \u2265 %d%%", c(25, 50, 75))),
          gate_lab = factor(fifelse(gate == "expected", GATES[1], GATES[2]), levels = GATES),
          stage_lab = fifelse(stage == "before", STAGES[1], STAGES[2]))]
sp[, row_lab := factor(paste0(gate_lab, "\n", stage_lab), levels = c(outer(GATES, STAGES, paste, sep = "\n")))]
fill_types <- scale_fill_manual(values = c(positive = unname(cmap_varType["AG"]), other = "grey90"), guide = "none")
lab_n <- function(n) format(n, big.mark = ",", trim = TRUE)
slug <- function(s) gsub("[^A-Za-z0-9]+", "_", tolower(s))
save_both <- function(p, stem, w, h) {
  save_figure_sd(p, stem, w, h, parts = "spectrum", grid = "y", script = SCRIPT, src_stem = STEM)
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}
spectrum_bars <- function(x, tot, n_size = 3, pct_size = 3, pct_face = "bold") {
  ggplot(x, aes(type, pct, fill = fifelse(type %in% c("AG", "TC"), "positive", "other"))) +
    geom_col(width = 0.75, colour = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 100 / 12, linetype = "dashed", colour = "grey50") +
    geom_text(data = x[type %in% c("AG", "TC")], aes(label = sprintf("%.1f", pct)),
              vjust = -0.4, size = pct_size, fontface = pct_face) +
    geom_text(data = tot, aes(x = 12.5, y = 98, label = lab), inherit.aes = FALSE, hjust = 1, vjust = 1,
              size = n_size) +
    fill_types + scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    theme_manuscript(base_size = BS, facets = TRUE) + grid_y
}

for (ds in unique(sp$dataset)) {
  for (tr in unique(sp[dataset == ds][order(truth_order), truth])) {
    x <- sp[dataset == ds & truth == tr]
    tot <- x[, .(n = sum(N), w1 = w1[1]), by = .(row_lab, t_lab, stage)]
    tot[, lab := fifelse(stage == "after" & !is.na(w1), sprintf("n = %s  (w1 = %d nt)", lab_n(n), w1),
                         sprintf("n = %s", lab_n(n)))]
    p <- spectrum_bars(x, tot) + facet_grid(row_lab ~ t_lab) +
      labs(x = "Mismatch type (read orientation)", y = "KO-depleted gate sites (%)",
           title = sprintf("%s ground truth from %s reads: mismatch spectrum", ds, tr),
           subtitle = "Dashed line: 12 types equally frequent. v9 positives = AG and TC after filtering.") +
      theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1))
    save_both(p, sprintf("GTDEF_bar_truth_spectrum12_%s_%s", slug(ds), slug(tr)), 12, 11)
  }
  y <- sp[dataset == ds & stage == "after" & t == "t50"]
  y[, truth := factor(truth, levels = unique(sp[order(truth_order), truth]))]
  tot <- y[, .(lab = paste0("n = ", lab_n(sum(N)))), by = .(truth, gate_lab)]
  p <- spectrum_bars(y, tot) + facet_grid(gate_lab ~ truth, labeller = label_wrap_gen(18)) +
    labs(x = "Mismatch type (read orientation)", y = "KO-depleted gate sites (%)",
         title = sprintf("%s ground truths by read set: mismatch spectrum after filtering", ds),
         subtitle = "Editing reduction in KO \u2265 50%. Dashed line: 12 types equally frequent.") +
    theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1))
  save_both(p, sprintf("GTDEF_bar_truth_spectrum12_%s_truths", slug(ds)), 16, 6.5)
}

# Response figure: headline truth of both datasets, t = 50%; rows = dataset, columns = gate x before / after filtering
r <- sp[truth_order == 1 & t == "t50"]
r[, `:=`(ds_lab = factor(c(Teratoma = "Human teratoma", Aorta = "Mouse aorta")[dataset],
                         levels = c("Human teratoma", "Mouse aorta")),
         col_lab = factor(fifelse(gate == "expected", "Power-gated", "Depth-gated"), levels = c("Power-gated", "Depth-gated")),
         stage_lab = factor(stage_lab, levels = STAGES))]
stopifnot(!anyNA(r$ds_lab))
tot <- r[, .(lab = sprintf("n = %s", lab_n(sum(N)))), by = .(ds_lab, col_lab, stage_lab)]
p <- spectrum_bars(r, tot, n_size = 4.2, pct_size = 4, pct_face = "plain") +
  ggh4x::facet_nested(ds_lab ~ col_lab + stage_lab) +
  labs(x = "Mismatch type (read orientation)",
       y = "Sites with mismatch rates\n\u2265 50% reduced in KO (%)") +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
        panel.spacing.y = unit(1, "lines"))
save_both(p, "GTDEF_bar_truth_spectrum12_response", 13, 6)
message("Done")
