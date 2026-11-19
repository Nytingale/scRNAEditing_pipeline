#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_config_factor_plot.R   (2026-10-15; draws from source data only since 2026-10-25)
# Response figure 13: paired F1 difference between two configurations that differ in one choice only, per
# tier and truth set (teratoma 6 KO + 3 REDIportal, held-out aorta 6 KO; v9, 76 configurations); one column per
# contrast (first column: each caller minus BCFtools), truth sets on the y-axis, points pool the tiers.
# Input : source_data/pm_gt/TERATOMA_v8_config_factor_delta_f1.pairs.tsv (TERATOMA_v8_config_factor_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_config_factor_delta_f1.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_config_factor_plot.R"
STEM <- "TERATOMA_v8_config_factor_delta_f1"

d <- read_part(STEM, "pairs")
y_max <- max(abs(d$delta)) * 1.05
# portrait layout: one row per contrast, one column per truth family, truth thresholds on the x-axis
truth_short <- function(x) {
  x <- sub(" reduction$", "", x)
  x <- sub("^≥ 3 variant reads,\n≥ 5% editing level$", "≥ 3 reads,\n≥ 5% level", x)
  sub("^≥ ([13]) variant (reads?)$", "≥ \\1 \\2", x)
}
d[, Truth := factor(truth_short(as.character(Truth)), rev(unique(truth_short(levels(Truth)))))]
fam_lab <- function(x) sub("^REDIportal-based", "REDIportal-\nbased", sub("\\(Human teratoma\\)", "(Human\nteratoma)", x))
d[, Family := factor(fam_lab(as.character(Family)), fam_lab(levels(Family)))]
choice_lab <- function(x) sub("^(GATK base) (quality)$", "\\1\n\\2", sub("^SNP removal$", "SNP\nremoval", x))
d[, Choice := factor(choice_lab(as.character(Choice)), choice_lab(levels(Choice)))]
d[, Detail := factor(gsub("\n", " ", as.character(Detail)), gsub("\n", " ", levels(Detail)))]
levels(d$Detail) <- sapply(levels(d$Detail), function(z) paste(strwrap(z, 16), collapse = "\n"))
BS <- 12

p <- ggplot(d, aes(Truth, delta)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(aes(colour = method), position = position_jitter(width = 0.22, height = 0, seed = 1), size = 0.9, alpha = 0.55) +
  geom_boxplot(fill = NA, colour = "black", width = 0.6, outlier.shape = NA, linewidth = 0.4) +
  ggh4x::facet_nested(Choice + Detail ~ Family, scales = "free_x", space = "free_x",
                      strip = ggh4x::strip_nested(text_y = ggh4x::elem_list_text(face = c("bold", "plain")),
                                                  by_layer_y = TRUE)) +
  scale_colour_manual_cvd(values = cmap_method, name = NULL,
                          guide = guide_legend(nrow = 1, override.aes = list(size = 3.5, alpha = 1))) +
  labs(x = NULL, y = "F1 difference (first minus second configuration)") +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-y_max, y_max)) +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_both +
  theme(legend.position = "bottom", legend.text = element_text(size = BS), strip.text.y = element_text(angle = 0, hjust = 0, size = BS - 1.5),
        strip.text.x = element_text(size = BS - 1), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = BS - 2, lineheight = 0.85),
        panel.spacing.y = unit(4, "pt"), panel.spacing.x = unit(5, "pt"),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey85"))
save_figure_sd(p, STEM, width = 9, height = 11.7, parts = c("pairs", "summary"), grid = "both", script = SCRIPT)
