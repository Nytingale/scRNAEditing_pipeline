#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_universe_celltype_robustness_plot.R   (split from GTDEF_gtdefinition_plot.R 2026-10-25; source data only)
# Response figure 2:
#   a  NTC mismatch spectrum at BQ/MQ 20 and 30        b  ADAR1 KO / NTC rate per mismatch type
#   c  NTC A-to-G allele fraction, in / not in dbSNP common
#   d  ADAR1 KO / NTC A-to-G fraction per NTC allele-fraction bin
#   e  truth sites per cell type (density over cell types, one curve per truth set)
#   f  Spearman's rho with the chosen-truth ranking, one threshold moved one step
#   g  rank of the chosen truth's top configuration
# Input : source_data/pm_gt/GTDEF_v8_universe_celltype_robustness.<part>.tsv
#         (GTDEF_v8_universe_celltype_robustness_source_data.R)
# Output: figures/pm_gt/GTDEF_v8_universe_celltype_robustness.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
suppressPackageStartupMessages({ library(dplyr); library(tibble) })
SCRIPT <- "publication/gt_definition/scripts/GTDEF_v8_universe_celltype_robustness_plot.R"
STEM <- "GTDEF_v8_universe_celltype_robustness"
PARTS <- c("a_spectrum", "b_ntc_vs_ko", "c_groups", "c_af_hist", "d_ko_response", "e_celltype", "f_rank_rho")

tn <- as_tibble(read_part(STEM, "a_spectrum"))
nk <- as_tibble(read_part(STEM, "b_ntc_vs_ko"))
grp_lab <- as_tibble(read_part(STEM, "c_groups"))
af <- as_tibble(read_part(STEM, "c_af_hist"))
kr <- as_tibble(read_part(STEM, "d_ko_response"))
dn <- read_part(STEM, "e_celltype")
sw <- read_part(STEM, "f_rank_rho")

GATE_COL <- c(expected = "#0072B2", depth = "#009E73")
gate_scale <- function(aes = "colour")
  scale_discrete_manual(aes, values = setNames(GATE_COL, V8_GATE_LAB), name = NULL)

# a-b: pileup quality
BS1 <- 13
p1a <- ggplot(tn, aes(type, pct_of_mismatches / 100, fill = key)) +
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
p1b <- ggplot(nk, aes(type, ko_over_ntc, colour = q)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual_cvd(values = cmap_bq_mq, name = NULL) +
  scale_y_continuous(limits = c(0.75, 1.05), breaks = seq(0.7, 1.1, 0.1)) +
  scale_x_discrete(labels = function(x) gsub(">", "", x)) +
  labs(x = "Mismatch type", y = "ADAR1 KO / NTC rate\n(per million bases)") +
  theme_manuscript(base_size = BS1) + grid_y + theme(legend.position = "top")

# c-d: site universe
BS2 <- 13
cols <- setNames(unname(cmap_dbsnp[grp_lab$group]), grp_lab$lab)
p2a <- ggplot(af, aes(100 * bin, 100 * frac, colour = g)) +
  geom_step(direction = "mid", linewidth = 0.9) +
  scale_colour_manual_cvd(values = cols, name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = "NTC A-to-G allele fraction (%)", y = "Share of sites (%)") +
  guides(colour = guide_legend(ncol = 1)) +
  theme_manuscript(base_size = BS2) + grid_y + theme(legend.position = "top")
p2b <- ggplot(kr, aes(bl, ko_over_ntc, colour = g, group = g)) +
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

# e: truth sites per cell type, nine-truth palette
NINE_COL <- setNames(unname(V8_TRUTH_COLOURS[V8_TRUTH_IDS]), levels(dn$Truth))
xb <- c(0, 10, V8_CT_MIN_SITES, 100, 1000, 10000, 1e5)
p_dn <- ggplot(dn, aes(y, colour = Truth, fill = Truth)) +
  geom_density(alpha = 0.06, linewidth = 0.9, bounds = c(0, Inf)) +
  geom_vline(xintercept = log10(V8_CT_MIN_SITES + 1), linetype = 2, linewidth = 0.6, colour = "black") +
  scale_x_continuous(breaks = log10(xb + 1), labels = format(xb, big.mark = ",", scientific = FALSE, trim = TRUE),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = NINE_COL, name = NULL) +
  scale_fill_manual(values = NINE_COL, name = NULL) +
  guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
  labs(x = "Number of truth sites in cell type", y = "Density of cell types") +
  theme_manuscript(base_size = 13) + grid_both +
  theme(legend.position = "inside", legend.position.inside = c(0.995, 0.995), legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key.size = unit(0.75, "lines"), legend.key.spacing.y = unit(0, "lines"),
        legend.margin = margin(2, 4, 2, 4))

# f-g: ranking under one-step threshold moves
sw_row <- function(p, ylab, x_text)
  p + facet_grid(. ~ Axis, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 14)) +
    scale_x_discrete(labels = function(x) sub("^[a-z_]+ ", "", x)) + gate_scale() + labs(x = NULL, y = ylab) +
    theme_manuscript(base_size = 13, facets = TRUE) + grid_y + theme(panel.spacing.x = unit(0.9, "lines")) +
    (if (!x_text) theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else theme())
p_sw_f <- sw_row(ggplot(sw, aes(x, rho, colour = g, group = g)) +
                   geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
                   geom_point(data = sw[centre == TRUE], shape = 21, size = 3.8, stroke = 1.1, fill = NA, colour = "black") +
                   geom_text(data = sw[centre == FALSE], aes(label = sprintf("%.3f", rho),
                             vjust = fifelse(gate == "depth", 1.8, -0.8)), size = 5, show.legend = FALSE,
                             colour = c(expected = "#003F62", depth = "#00573F")[sw[centre == FALSE]$gate]) +
                   scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)),
                 "Spearman's \u03C1 with the\nchosen-truth ranking", x_text = FALSE) +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1))
p_sw_g <- sw_row(ggplot(sw, aes(x, centre_top_rank, colour = g, group = g)) +
                   geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
                   geom_point(data = sw[centre == TRUE], shape = 21, size = 3.8, stroke = 1.1, fill = NA, colour = "black") +
                   scale_y_reverse(limits = c(max(3, max(sw$centre_top_rank)) + 0.3, 0.7),
                                   breaks = 1:max(3, max(sw$centre_top_rank))),
                 "Rank of the chosen\ntruth's top configuration", x_text = TRUE) +
  theme(legend.position = "none", strip.text.x = element_blank(), strip.background = element_blank())

inside <- function(x, y, jx = x, jy = y) theme(legend.position = "inside", legend.position.inside = c(x, y),
                                              legend.justification = c(jx, jy), legend.background = element_blank())
# text sizes per panel (a figure-wide & theme() would re-enable the blank x labels of f)
UE_TXT <- 17
big <- function(p, x_text = TRUE)
  p + theme(axis.text.x = if (x_text) element_text(size = UE_TXT) else element_blank(),
            axis.text.y = element_text(size = UE_TXT), axis.title.x = element_text(size = UE_TXT + 1),
            axis.title.y = element_text(size = UE_TXT + 1), strip.text = element_text(size = UE_TXT),
            legend.text = element_text(size = UE_TXT - 1), legend.title = element_text(size = UE_TXT - 1))
r1 <- patchwork::wrap_plots(big(p1a) + theme(axis.text.x = element_text(size = UE_TXT - 3)), big(p2a + inside(0.5, 0.98)),
                            big(p2b + guides(colour = guide_legend(ncol = 1)) + inside(0.98, 0.03) +
                                  theme(axis.text.x = element_text(angle = 35, hjust = 1))),
                            nrow = 1, widths = c(1.3, 1, 1.1))
r2 <- patchwork::wrap_plots(patchwork::free(big(p1b + inside(0.98, 0.03)) +
                                              theme(axis.text.x = element_text(size = UE_TXT - 3, angle = 90, vjust = 0.5, hjust = 1))),
                            big(p_dn) + theme(legend.text = element_text(size = UE_TXT - 4)),
                            nrow = 1, widths = c(1.3, 2.1))
fig_ue <- patchwork::wrap_plots(r1, r2, big(p_sw_f, x_text = FALSE) + labs(y = "Spearman's \u03C1 with\nchosen ranking"),
                                big(p_sw_g) + labs(y = "Rank of chosen\ntop configuration") +
                                  theme(axis.text.x = element_text(size = UE_TXT - 2)), ncol = 1,
                                heights = c(1, 1.2, 0.75, 0.55)) +
  patchwork::plot_annotation(tag_levels = list(c("a", "c", "d", "b", "e", "f", "g"))) &
  theme(plot.tag = element_text(face = "bold", size = 26))
save_figure_sd(fig_ue, STEM, width = 16, height = 21, parts = PARTS, grid = "y", script = SCRIPT)
