#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_informative_ko_plot.R   (split from GTDEF_gtdefinition_plot.R 2026-10-25; source data only)
# Response figure 1: the informative KO arm.
#   a  e^-m, the chance of a false complete loss, at m = 3 (chosen) and 4.6 (theory)
#   b  sites with 0 variant reads in the KO vs m: non-AG/TC sites follow e^-m, A-to-G sites lie above it
#   c  KO coverage each gate requires vs NTC editing level (power-gated: 3 / p; depth-gated: 10) (theory)
#   d  power-gated: sites vs m          e  depth-gated: sites vs KO coverage
# Printed 6.5 in wide in the response: larger text than the supplementary figure.
# Input : source_data/pm_gt/GTDEF_v8_informative_ko.{b_null_fit,d_m_curve,e_cov_curve}.tsv
#         (GTDEF_v8_informative_ko_source_data.R)
# Output: figures/pm_gt/GTDEF_v8_informative_ko.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble) })
source(file.path(SCRIPT_DIR, "GTDEF_v8_m_null_panels.R"))
SCRIPT <- "publication/gt_definition/scripts/GTDEF_v8_informative_ko_plot.R"
STEM <- "GTDEF_v8_informative_ko"

zf <- as_tibble(read_part(STEM, "b_null_fit"))
mcv <- as_tibble(read_part(STEM, "d_m_curve"))
kcv <- as_tibble(read_part(STEM, "e_cov_curve"))

BSP <- 17; LS <- 18
GATE_COL <- c(expected = "#0072B2", depth = "#009E73")

# d-e: sites as the informative-KO threshold of each gate moves
cs_lab0 <- c(n_true_t0 = "\u2265 0% editing reduction in KO", n_true_t25 = "\u2265 25% editing reduction in KO",
             n_true = "\u2265 50% editing reduction in KO", n_true_t75 = "\u2265 75% editing reduction in KO")
cs_lt <- c("solid", "solid", "solid", "solid")
curve_panel <- function(d, x, chosen, xlab, ttl, x_step, cols, sl) {
  keep = intersect(names(sl), names(d))
  dl = d %>% pivot_longer(all_of(keep), names_to = "s", values_to = "n") %>%
    mutate(s = factor(sl[s], levels = sl[keep]))
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
CS_EXP <- c("black", "#99C7E0", "#0072B2", "#003F62"); CS_DEP <- c("black", "#99D8C7", "#009E73", "#00573F")
leg_bl <- theme(legend.position = "inside", legend.position.inside = c(0.01, 0.01), legend.justification = c(0, 0),
                legend.background = element_blank(), legend.key.size = unit(1, "lines"),
                legend.text = element_text(size = BSP - 2))
leg_room <- function(d) expand_limits(y = log10(min(d$n_true_t75)) - 0.55)
p_h <- curve_panel(mcv, "ko_min_expected", 3, "m = expected edited KO reads", V8_GATE_LAB2[["expected"]], 1,
                   cols = CS_EXP, sl = cs_lab0) + leg_bl + leg_room(mcv)
p_i <- curve_panel(kcv, "ko_min_cov", 10, "KO coverage (reads)", V8_GATE_LAB2[["depth"]], 5,
                   cols = CS_DEP, sl = cs_lab0) + leg_bl + leg_room(kcv)

# a: exp(-m), the chance of a false complete loss, marked at m = 3 (chosen) and 4.6
M_PTS <- c(3, 4.6)
mpt <- tibble(m = M_PTS, p0 = exp(-M_PTS)) %>%
  mutate(chosen = abs(m - 3) < 1e-9, hj = if_else(m > 4, 1.12, -0.12), vj = if_else(m > 4, 1.4, -0.5),
         lab = sprintf("%s%% chance: m = %s", formatC(100 * p0, format = "fg", digits = 1), as.character(m)))
p_fa <- ggplot(tibble(m = seq(0, 8, 0.02)) %>% mutate(p0 = exp(-m)), aes(m, p0)) +
  geom_hline(yintercept = 0.05, linetype = 3, colour = "grey55") +
  geom_line(linewidth = 1.1) +
  geom_point(data = mpt, aes(colour = chosen), size = 3.2) +
  geom_text(data = mpt, aes(label = lab, colour = chosen, hjust = hj, vjust = vj), size = 5.6, lineheight = 0.9) +
  scale_colour_manual(values = c(`TRUE` = GATE_COL[["expected"]], `FALSE` = "grey25"), guide = "none") +
  scale_y_log10(breaks = c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001),
                labels = function(x) paste0(100 * x, "%")) +
  scale_x_continuous(breaks = 0:8) +
  coord_cartesian(xlim = c(0, 8), ylim = c(1e-3, 1)) +
  labs(x = "m = expected edited KO reads\n(KO coverage \u00D7 NTC editing level)",
       y = "Chance of 0 edited KO reads\nwhen editing is unchanged (e\u207B\u1D50)", title = "False complete loss") +
  theme_manuscript(base_size = BSP) + grid_both
# c: KO coverage each gate requires vs NTC editing level (m = 3: 3 / p; depth: 10), from the 5% NTC editing floor
F2_COL <- GATE_COL
kn <- tibble(p = 10^seq(log10(0.05), 0, length.out = 200)) %>% mutate(expected = 3 / p, depth = 10) %>%
  pivot_longer(c(expected, depth), names_to = "gate", values_to = "n")
kpt <- tibble(p = c(0.05, 0.1, 0.5, 1), n = 3 / c(0.05, 0.1, 0.5, 1)) %>%
  mutate(lab = sprintf("%g%%: %g reads", 100 * p, n), hj = if_else(p >= 0.5, 1.1, -0.12), vj = if_else(p >= 0.5, 1.6, -0.5))
p_kn <- ggplot(kn, aes(100 * p, n, colour = gate)) +
  geom_line(linewidth = 1.1, show.legend = FALSE) +
  geom_point(data = kpt, aes(100 * p, n), inherit.aes = FALSE, size = 2.8, colour = F2_COL[["expected"]]) +
  geom_text(data = kpt, aes(100 * p, n, label = lab, hjust = hj, vjust = vj), inherit.aes = FALSE, size = 5.4,
            colour = F2_COL[["expected"]]) +
  annotate("text", x = 106, y = 11.5, label = "Depth-gated: \u2265 10 KO reads", hjust = 1, vjust = 0, size = 5.4,
           colour = F2_COL[["depth"]]) +
  annotate("text", x = 35, y = 45, label = "Power-gated: KO coverage\n\u2265 3 / NTC editing level",
           hjust = 0, vjust = 0, size = 5.4, colour = F2_COL[["expected"]], lineheight = 0.9) +
  scale_colour_manual(values = F2_COL, guide = "none") +
  scale_x_continuous(breaks = c(5, 25, 50, 75, 100), limits = c(0, 108)) +
  scale_y_continuous(breaks = c(3, 10, 20, 30, 40, 50, 60), limits = c(0, 66)) +
  labs(x = "NTC editing level (%)", y = "KO coverage required (reads)", title = "KO coverage required") +
  theme_manuscript(base_size = BSP) + grid_both

ttl <- theme(plot.title = element_text(size = BSP + 2, face = "bold", hjust = 0.5))
p_null <- m_null_zero_panel(zf, BSP) +
  theme(legend.position = "bottom", legend.text = element_text(size = BSP + 1)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, keywidth = unit(2.2, "lines")))
row1 <- patchwork::wrap_plots(p_fa, p_null, nrow = 1, widths = c(1, 2.3)) & ttl
row2 <- patchwork::wrap_plots(p_kn, p_h, p_i, nrow = 1) & ttl
r_f1 <- patchwork::wrap_plots(row1, row2, ncol = 1, heights = c(1, 1)) +
  patchwork::plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = LS))
save_figure_sd(r_f1, STEM, width = 17, height = 12, parts = c("b_null_fit", "d_m_curve", "e_cov_curve"), grid = "y",
               script = SCRIPT)
