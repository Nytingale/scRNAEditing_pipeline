#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_spectrum12_plot.R
#
# Mismatch spectrum of the v6 ground truth: the frozen v6 rules applied to all
# 12 mismatch types (RNA strand), for every depletion threshold d of the step-6
# sweep (25-90% in 5% steps).
#   GTDEF_bar_truth_spectrum12_d<NN>   one figure per d: ground-truth sites per type
#   GTDEF_bar_truth_spectrum12_all_d   all d, % of ground-truth sites per type
#   GTDEF_line_truth_spectrum12        A>G vs the other 11 types across d
#
# INPUT : human/pm/gt_eval_v6/spectrum12/eval/spectrum12_by_d.tsv, parity_AG_d050.tsv
#         (human/pm/scripts/run_gt_spectrum12_v6.sh)
# OUTPUT: publication/gt_definition/{figures,source_data}/<DATASET>/
# Usage : Rscript GTDEF_spectrum12_plot.R [DATASET]
# ============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "GTDEF_common.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
argv = commandArgs(trailingOnly = TRUE)
gtdef_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT = caller_script()

EV = file.path(GT_V6_DIR, "spectrum12", "eval")
SRC = file.path(EV, "spectrum12_by_d.tsv")
PAR = file.path(EV, "parity_AG_d050.tsv")
BS = 13
TYPES = all_vartypes()
sp = fread(SRC)
sp = merge(CJ(type = TYPES, d = unique(sp$d)), sp, by = c("type", "d"), all.x = TRUE)
for (cn in c("n_candidates", "n_gate", "n_true", "n_not_depleted")) sp[is.na(get(cn)), (cn) := 0L]
sp[, `:=`(type = factor(type, levels = TYPES), pct = 100 * n_true / sum(n_true)), by = d]
sp[, d_pct := as.integer(round(100 * d))]
tot = sp[, .(n_total = sum(n_true), n_ag = sum(n_true[type == "AG"])), by = d_pct]
tot[, pct_ag := 100 * n_ag / n_total]
par = if (file.exists(PAR)) fread(PAR) else NULL
if (!is.null(par)) message("parity AG d=0.50: v6 ", par$n_v6, ", spectrum12 ", par$n_spectrum12,
                           ", shared ", par$n_shared)

fill_types = scale_fill_manual(values = cmap_varType, guide = "none")
x_types = scale_x_discrete(drop = FALSE)

# ---- one figure per d ------------------------------------------------------
for (dp in sort(unique(sp$d_pct))) {
  x = sp[d_pct == dp]
  p = ggplot(x, aes(type, n_true, fill = type)) +
    geom_col(width = 0.75, colour = "grey20", linewidth = 0.2) +
    geom_text(aes(label = format(n_true, big.mark = ",", trim = TRUE)), vjust = -0.4, size = 3.4) +
    fill_types + x_types +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = shorten_axis) +
    labs(x = "Mismatch type", y = "Ground-truth sites",
         title = sprintf("Editing reduction in KO \u2265 %d%%", dp)) +
    theme_manuscript(base_size = BS) + grid_y
  save_figure(p, sprintf("GTDEF_bar_truth_spectrum12_d%02d", dp), 7.5, 4.8,
              sources = c(SRC, PAR), data = x, grid = "y", script = SCRIPT)
}

# ---- all d, % per type -----------------------------------------------------
lab_d = setNames(sprintf("\u2265 %d%%  (n = %s)", tot$d_pct, format(tot$n_ag, big.mark = ",", trim = TRUE)), tot$d_pct)
sp[, d_lab := factor(lab_d[as.character(d_pct)], levels = lab_d[as.character(sort(tot$d_pct))])]
p_all = ggplot(sp, aes(type, pct, fill = type)) +
  geom_col(width = 0.75, colour = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = 100 / 12, linetype = "dashed", colour = "grey50") +
  geom_text(data = sp[type == "AG"], aes(label = sprintf("%.1f", pct)),
            vjust = -0.5, size = 3.4, fontface = "bold", colour = "black") +
  facet_wrap(~ d_lab, ncol = 4) +
  fill_types + x_types +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Mismatch type", y = "Ground-truth sites (%)",
       title = "Mismatch types among ground-truth sites, by editing reduction in KO") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y +
  theme(axis.text.x = element_text(size = 9))
save_figure(p_all, "GTDEF_bar_truth_spectrum12_all_d", 16, 12,
            sources = c(SRC, PAR), data = sp, grid = "y", script = SCRIPT)

# ---- A>G vs the other 11 types across d -----------------------------------
ln = melt(tot, id.vars = "d_pct", measure.vars = c("n_ag", "n_total"), variable.name = "set", value.name = "n")
ln[set == "n_total", n := n - tot$n_ag[match(d_pct, tot$d_pct)]]
ln[, set := factor(fifelse(set == "n_ag", "AG", "Other 11 types"), levels = c("AG", "Other 11 types"))]
cmap_ag = c("AG" = cmap_varType[["AG"]], "Other 11 types" = "grey50")
p_n = ggplot(ln, aes(d_pct, n, colour = set)) +
  geom_vline(xintercept = 100 * GT_DEPLETION, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = cmap_ag, name = NULL) +
  scale_x_continuous(breaks = seq(25, 90, 10)) +
  scale_y_continuous(labels = shorten_axis) +
  labs(x = "Editing reduction in KO (%)", y = "Ground-truth sites") +
  theme_manuscript(base_size = BS) + grid_both + theme(legend.position = "top")
p_pct = ggplot(tot, aes(d_pct, pct_ag)) +
  geom_vline(xintercept = 100 * GT_DEPLETION, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8, colour = cmap_varType[["AG"]]) +
  geom_point(size = 2, colour = cmap_varType[["AG"]]) +
  scale_x_continuous(breaks = seq(25, 90, 10)) +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(x = "Editing reduction in KO (%)", y = "AG among ground-truth sites (%)") +
  theme_manuscript(base_size = BS) + grid_both
p_line = cowplot::plot_grid(p_n, p_pct, labels = c("a", "b"), label_size = 16, nrow = 1, align = "h", axis = "b")
save_figure(p_line, "GTDEF_line_truth_spectrum12", 11, 4.8,
            sources = c(SRC, PAR), data = merge(ln, tot, by = "d_pct"), grid = "both", script = SCRIPT)
message("Done")
