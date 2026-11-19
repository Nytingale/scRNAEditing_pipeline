#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_caller_margin_plot.R   (2026-11-02)
# Margin of the winning configuration over the best configuration of a different caller, per truth set
# (mean F1 across tiers), absolute difference only (2026-11-19).
# Input : source_data/pm_gt/TERATOMA_v8_caller_margin.margin.tsv (TERATOMA_v8_caller_margin_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_caller_margin.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_caller_margin_plot.R"
STEM <- "TERATOMA_v8_caller_margin"

m <- read_part(STEM, "margin")
MEAS <- c(abs = "Difference in mean F1", pct = "Difference as % of the winner's mean F1")
d <- m[, .(Truth, Family, Winner, Other, measure = MEAS[["abs"]], v = margin, lab = sprintf("%+.3f", margin))]
d[, `:=`(measure = factor(measure, MEAS),
         lab = paste0(lab, "  ", Winner, " vs ", Other))]
# labels start right of the bar (or of 0 for negative bars); room for them on the right, per column
lim <- d[, .(hi = max(v), lo = min(min(v), 0)), by = measure][, span := hi - lo]
pad <- rbind(lim[, .(measure, v = hi + 1.15 * span)], lim[, .(measure, v = 1.3 * lo)])
ft <- unique(d[, .(Family, Truth)])[, .SD[1], by = Family]
pad <- rbindlist(lapply(seq_len(nrow(ft)), function(i) cbind(pad, ft[i])))

assert_cvd_safe(cmap_method, "scale_fill_manual")
p <- ggplot(d, aes(v, Truth)) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_col(aes(fill = Winner), colour = "grey15", linewidth = 0.3, width = 0.7) +
  geom_text(aes(x = pmax(v, 0), label = lab), hjust = -0.05, size = 3.6) +
  geom_blank(data = pad, aes(v, Truth)) +
  facet_grid(Family ~ ., scales = "free", space = "free_y") +
  scale_fill_manual(values = cmap_method, name = "Caller ranked first") +
  labs(x = "\u0394F1: winning configuration minus the best\nconfiguration of a different caller (mean F1)", y = NULL) +
  theme_manuscript(base_size = 13, facets = TRUE) + grid_x +
  theme(legend.position = "bottom")
save_figure_sd(p, STEM, width = 8.5, height = 5.5, parts = "margin", grid = "x", script = SCRIPT)
