#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_negatives_plot.R   (2026-10-15; draws from source data only since 2026-10-25)
# Response figure 3: pooled NTC editing level of the informative negatives vs the positives (t = 50%), per
# coverage gate, and the share of each class in REDIportal.
# Input : source_data/pm_gt/TERATOMA_v8_negatives_ntc_editing.{sites,labels}.tsv
#         (TERATOMA_v8_negatives_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_negatives_ntc_editing.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_negatives_plot.R"
STEM <- "TERATOMA_v8_negatives_ntc_editing"

d <- read_part(STEM, "sites")
lab <- read_part(STEM, "labels")
cmap_class <- setNames(c("#2171B5", "#CB181D"), levels(d$Class))
BIN <- 0.05
d[, bin := pmin(floor(edit / BIN), 1 / BIN - 1)]
y_top <- d[, .N, by = .(Gate, Class, bin)][, pct := 100 * N / sum(N), by = .(Gate, Class)][, max(pct)]
lab[, y := y_top * c(1.22, 1.12)[as.integer(Class)]]

p <- ggplot(d, aes(edit, fill = Class)) +
  geom_histogram(aes(y = after_stat(100 * density * BIN)), binwidth = BIN, boundary = 0, position = "identity",
                 alpha = 0.55, colour = NA) +
  geom_text(data = lab, aes(x = 0.5, y = y, label = txt, colour = Class), hjust = 0.5, size = 4.2, inherit.aes = FALSE,
            show.legend = FALSE) +
  facet_wrap(~ Gate, nrow = 1) +
  scale_fill_manual(values = cmap_class, name = NULL) +
  scale_colour_manual(values = cmap_class, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, y_top * 1.28), expand = expansion(mult = c(0, 0))) +
  labs(x = "Editing level in merged NTC arm", y = "Percentage of A-to-G sites (%)") +
  theme_manuscript(base_size = 14, facets = TRUE) + grid_x
save_figure_sd(p, STEM, width = 12, height = 5.5, parts = c("sites", "labels"), grid = "x", script = SCRIPT)
