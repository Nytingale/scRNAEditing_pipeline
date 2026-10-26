#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_sweep_rank_plot.R   (2026-10-14; layout as GTDEF_supp_gt_thresholds_v8 panel e, 2026-10-15)
# One truth threshold at a time moved one step down and up from the chosen value (circled), both coverage gates,
# t = 50% KO editing reduction (as GTDEF panel e; t = 25% / 75% in the source data):
#   top     Spearman rho of the overall configuration ranking vs the ranking under the chosen truth
#   bottom  rank, under the neighbouring truth, of the configuration ranked first under the chosen truth
# Gate-specific thresholds (expected KO reads: expected gate; KO coverage: depth gate) only for their gate.
# Input : human/pm/gt_eval_v8/sweep_rank/sweep_rank_compare<V8_SFX>.tsv (human/scripts/gt_v8_sweep_rank_compare.R;
#         headline _pwcov, plotting/gt_v8_paths.R)
# Output: figures/pm_gt/TERATOMA_v8_sweep_rank_rho.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_sweep_rank_plot.R"

SRC <- file.path(V8_DIR, "sweep_rank", paste0("sweep_rank_compare", V8_SFX, ".tsv"))
stopifnot(file.exists(SRC))
# chosen values: gt_eval_v6/FROZEN_RULES.tsv (human/scripts/gt_robustness_v8.R CENTRE)
CENTRE <- c(min_cov = "10", min_edit_reads = "3", min_ctrl_edit = "0.05", min_cov_teratoma = "10",
            rule = "evaluable", ko_min_expected = "3", ko_min_cov = "10")
AXIS_LAB <- c(min_cov = "NTC coverage", min_edit_reads = "NTC edited reads", min_ctrl_edit = "NTC editing level",
              min_cov_teratoma = "Per-sample NTC coverage", rule = "Reproducibility across samples",
              ko_min_expected = "Expected edited KO reads (m)", ko_min_cov = "KO coverage")
GATE_ONLY <- c(ko_min_expected = "expected", ko_min_cov = "depth")
GATE_COL <- setNames(c("#0072B2", "#009E73"), V8_GATE_LAB[c("expected", "depth")])

all_t <- fread(SRC)[truth != "centre"]
all_t[, value := substring(truth, nchar(parameter) + 2L)]
stopifnot(all(startsWith(all_t$truth, paste0(all_t$parameter, "_"))))
stopifnot(all(all_t$parameter %in% names(CENTRE)))
all_t <- all_t[!(parameter %in% names(GATE_ONLY)) | gate == GATE_ONLY[parameter]]
cen <- unique(all_t[, .(gate, t, parameter)])[, `:=`(value = CENTRE[parameter], rho = 1, centre_top_rank = 1L)]
d <- rbind(all_t[, .(gate, t, parameter, value, rho, centre_top_rank)], cen)[, centre := value == CENTRE[parameter]]
print(d[centre == FALSE, .(min_rho = round(min(rho), 3), max_rank_of_top = max(centre_top_rank)), by = t])

d <- d[t == "t50"]
num <- suppressWarnings(as.numeric(d$value))
RULE_ORDER <- c(strict2 = 1, evaluable = 2, strict3 = 3)
d[, ord := fifelse(parameter == "rule", RULE_ORDER[value], num)]
d[, lab := fcase(parameter == "min_ctrl_edit", paste0(round(100 * num), "%"),
                 parameter == "rule" & value == "strict2", "\u22652\nof 4",
                 parameter == "rule" & value == "strict3", "\u22653\nof 4",
                 parameter == "rule", "all\nbut 1", default = value)]
setorder(d, parameter, ord)
d[, x := factor(paste(parameter, lab), unique(paste(parameter, lab)))]
d[, `:=`(Axis = factor(AXIS_LAB[parameter], AXIS_LAB), Gate = factor(V8_GATE_LAB[gate], V8_GATE_LAB[c("expected", "depth")]))]
relab <- function(x) sub("^[a-z_]+ ", "", x)

row <- function(p, ylab, x_text)
  p + facet_grid(. ~ Axis, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 14)) +
    scale_x_discrete(labels = relab) + scale_colour_manual(values = GATE_COL, name = NULL) +
    labs(x = NULL, y = ylab) + theme_manuscript(base_size = 14, facets = TRUE) + grid_y +
    theme(panel.spacing.x = unit(0.9, "lines")) +
    (if (!x_text) theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else theme())
p_top <- row(ggplot(d, aes(x, rho, colour = Gate, group = Gate)) +
               geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
               geom_point(data = d[centre == TRUE], shape = 21, size = 3.8, stroke = 1.1, fill = NA, colour = "black") +
               geom_text(data = d[centre == FALSE], aes(label = sprintf("%.3f", rho), vjust = fifelse(gate == "depth", 1.8, -0.8)),
                         size = 4.3, show.legend = FALSE) +
               scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)),
             "Spearman's \u03C1 with the\nchosen-truth ranking", x_text = FALSE) +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1))
p_bot <- row(ggplot(d, aes(x, centre_top_rank, colour = Gate, group = Gate)) +
               geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
               geom_point(data = d[centre == TRUE], shape = 21, size = 3.8, stroke = 1.1, fill = NA, colour = "black") +
               scale_y_reverse(limits = c(max(3, max(d$centre_top_rank)) + 0.3, 0.7), breaks = 1:max(3, max(d$centre_top_rank))),
             "Rank of the chosen\ntruth's top configuration", x_text = TRUE) +
  theme(legend.position = "none", strip.text.x = element_blank(), strip.background = element_blank())
p <- (p_top / p_bot) + plot_layout(heights = c(1.2, 1))
save_figure(p, "TERATOMA_v8_sweep_rank_rho", width = 15, height = 7, sources = SRC,
            data = d[, .(gate, t, parameter, value, centre, rho, centre_top_rank)], grid = "y", script = SCRIPT)
