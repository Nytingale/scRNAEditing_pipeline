#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_metric_case_plot.R   (2026-10-21; draws from source data only since 2026-10-25)
# Response figure 9: truth-free measures of the 58 configurations, standard downstream filters.
#   a  rank by number of A-to-G sites vs rank by A-to-G percentage; outlined = rank shift >= 20 (dotted)
#   b  mean recall and mean precision vs mean F1
#   c  Spearman's rho of each measure with mean F1 and with the KO response of the whole call set
# Input : source_data/pm_gt/TERATOMA_v8_metric_case.{a_ranks,c_rho}.tsv (TERATOMA_v8_metric_case_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_metric_case.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_metric_case_plot.R"
STEM <- "TERATOMA_v8_metric_case"
SHIFT <- 20

s0 <- read_part(STEM, "a_ranks")
rc <- read_part(STEM, "c_rho")
shape_proc <- setNames(c(21, 24, 22, 23, 25)[seq_along(levels(s0$proc))], levels(s0$proc))
rho <- function(x, y) cor(x, y, method = "spearman")

BASE <- 17
F1_LAB <- "Mean F1 across 3 tiers\nand 6 KO truth sets"
cor_lab <- function(x, y) sprintf("Pearson's r = %.2f\n%s", cor(x, y), spearman_lab(rho(x, y)))
look <- list(theme_manuscript(base_size = BASE, facets = TRUE), grid_both, theme(aspect.ratio = 1))
scales_pts <- list(
  scale_fill_manual(values = cmap_method, name = NULL,
                    guide = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, size = 4.5, alpha = 0.9))),
  scale_shape_manual(values = shape_proc, name = NULL,
                     guide = guide_legend(nrow = 2, order = 2, override.aes = list(fill = "grey60", size = 4.5, alpha = 0.9))))
yup <- scale_y_continuous(limits = c(0, NA), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.05, 0.3)))

BR <- c(1, seq(10, 50, 10), 58)
p_top_a <- ggplot(copy(s0)[, panel := "Main analysis"], aes(r_ag, r_pct)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey20", linewidth = 0.6) +
  geom_abline(slope = 1, intercept = c(-SHIFT, SHIFT), linetype = "dotted", colour = "grey15", linewidth = 0.9) +
  geom_point(data = copy(s0[abs(shift) < SHIFT])[, panel := "Main analysis"], aes(fill = method, shape = proc),
             colour = "grey40", stroke = 0.3, size = 2.6, alpha = 0.3) +
  geom_point(data = copy(s0[abs(shift) >= SHIFT])[, panel := "Main analysis"], aes(fill = method, shape = proc),
             colour = "black", stroke = 0.9, size = 3.4, alpha = 0.95) +
  annotate("text", x = 58, y = -Inf, label = cor_lab(s0$r_ag, s0$r_pct), hjust = 0, vjust = 1.2, size = 5.2,
           lineheight = 0.9) +
  scale_x_reverse(breaks = BR, expand = expansion(add = 1)) +
  scale_y_reverse(breaks = BR, expand = expansion(add = c(1, 14))) +
  facet_wrap(~ panel) + scales_pts + look + guides(fill = "none", shape = "none") +
  labs(x = "Rank by number of\nA-to-G sites", y = "Rank by A-to-G\npercentage")

dl <- rbind(s0[, .(method, proc, mean_f1, metric = "Recall", value = recall)],
            s0[, .(method, proc, mean_f1, metric = "Precision", value = precision)])
dl[, metric := factor(metric, c("Recall", "Precision"))]
rl <- dl[, .(lab = cor_lab(value, mean_f1)), by = metric]
ll <- dl[, {fit <- lm(mean_f1 ~ value); xg <- seq(min(value), max(value), length.out = 100)
            .(x = xg, y = coef(fit)[1] + coef(fit)[2] * xg)}, by = metric]
p_top_b <- ggplot(dl, aes(value, mean_f1)) +
  geom_line(data = ll, aes(x, y), linetype = "dashed", colour = "grey20", linewidth = 0.7) +
  geom_point(aes(fill = method, shape = proc), colour = "grey30", stroke = 0.3, size = 2.8, alpha = 0.6) +
  geom_text(data = rl, aes(x = 0, y = Inf, label = lab), hjust = 0, vjust = 1.2, size = 5.2, lineheight = 0.9,
            inherit.aes = FALSE) +
  facet_wrap(~ metric, nrow = 1) + scales_pts + yup +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) + look +
  labs(x = "Mean across 3 tiers and 6 KO truth sets", y = F1_LAB)

p_top_c <- ggplot(rc, aes(measure, rho)) +
  geom_col(aes(fill = target), width = 0.7, colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", rho)), vjust = -0.4, size = 5) +
  facet_wrap(~ target, nrow = 1) +
  scale_fill_manual(values = setNames(c("grey65", "#4DAF4A"), levels(rc$target)), guide = "none") +
  scale_y_continuous(limits = c(0, 1.12), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0))) +
  theme_manuscript(base_size = BASE, facets = TRUE) + grid_both +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, lineheight = 0.85)) +
  labs(x = NULL, y = "Spearman's \u03c1")
top <- (p_top_a | p_top_b | p_top_c) + plot_layout(widths = c(1, 2, 1.7))
fig_theme <- theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(2, "pt"),
                   plot.tag = element_text(face = "bold", size = 20))
p <- top + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & fig_theme
save_figure_sd(p, STEM, width = 17, height = 7.4, parts = c("a_ranks", "c_rho"), grid = "both", script = SCRIPT)
