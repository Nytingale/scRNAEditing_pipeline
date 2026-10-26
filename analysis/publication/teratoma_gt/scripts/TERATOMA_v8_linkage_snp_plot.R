#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_linkage_snp_plot.R   (2026-10-17; draws from source data only since 2026-10-25)
# Review R3.9 (linked allele-specific markers): linking candidate editing sites to heterozygous PGP1 SNVs, NTC arm.
#   _allpanels  (a) sites with a heterozygous SNV within one read length; (b) site x unit pairs with >= 1 and >= 10
#               linked molecules in pseudobulks and single cells; (c) NTC cells needed for 10 linked molecules;
#               (d) site x unit pairs testable by bulk allele-specific editing criteria
#   _cells      response figure 18: KO-based positives with an edited linked molecule and with both SNV alleles covered
# Input : source_data/pm_gt/TERATOMA_v8_linkage_snp.{a_sites,b_pairs,c_cells,c_refs,d_ased,ko_steps}.tsv
#         (TERATOMA_v8_linkage_snp_source_data.R)
# Output: figures/pm_gt/TERATOMA_v8_linkage_snp_cells{,_allpanels}.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_plot_common.R"))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_linkage_snp_plot.R"
STEM <- "TERATOMA_v8_linkage_snp"

pa_d <- read_part(STEM, "a_sites")
pb_d <- read_part(STEM, "b_pairs")
pc_d <- read_part(STEM, "c_cells")
refs <- read_part(STEM, "c_refs")
pd_d <- read_part(STEM, "d_ased")
fn <- read_part(STEM, "ko_steps")

LEVELS <- levels(fn$Level)
cmap_set <- setNames(c("#2171B5", "grey45"), levels(pa_d$Set))
pa <- ggplot(pa_d, aes(what, pct, fill = Set)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%s of\n%s", format(n_hit, big.mark = ","), format(n, big.mark = ","))),
            position = position_dodge(0.8), vjust = -0.3, size = 3.4, lineheight = 0.9) +
  scale_fill_manual(values = cmap_set, guide = "none") +
  scale_y_continuous(limits = c(0, 20), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage of sites (%)") +
  theme_manuscript(base_size = 13, axis_text_size = 11.5) + grid_y

short_n <- label_number(scale_cut = cut_short_scale(), accuracy = 0.1, drop0trailing = TRUE)
x_tot <- dcast(unique(pb_d[, .(x, set, n_pairs)]), x ~ set, value.var = "n_pairs")
x_lab <- setNames(sprintf("%s\n%s | %s", x_tot$x, short_n(x_tot$ko), short_n(x_tot$cand)), x_tot$x)
pb <- ggplot(pb_d, aes(x, pct, colour = Set, shape = Threshold)) +
  geom_point(size = 3.2) +
  geom_line(data = pb_d[as.integer(x) <= 4], aes(group = interaction(Set, Threshold)), linewidth = 0.6) +
  geom_vline(xintercept = 3.5, colour = "grey60", linewidth = 0.4) +
  annotate("text", x = c(2, 5.5), y = 300, label = c("Pseudobulks", "Single NTC cells"), size = 4.2,
           colour = "grey30") +
  scale_x_discrete(labels = x_lab) +
  scale_colour_manual(values = cmap_set, name = NULL) +
  scale_shape_manual(values = c(16, 17), name = "Linked molecules per site and unit") +
  scale_y_log10(limits = c(1e-6, 500), breaks = 10^(-6:2),
                labels = function(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 6))) +
  labs(x = "Possible site \u00d7 unit pairs: KO-based positives | all candidate A-to-G sites",
       y = "Site \u00d7 unit pairs with linked molecules (%)") +
  theme_manuscript(base_size = 13, axis_text_size = 11) + grid_y

refs <- refs[, .(x, hj = c(1.05, -0.05), lab)]
pc <- ggplot(pc_d, aes(cells_needed, pct_sites, colour = Set)) +
  geom_vline(data = refs, aes(xintercept = x), linetype = "dashed", colour = "grey40") +
  geom_text(data = refs, aes(x = x, y = 14.6, label = lab, hjust = hj), inherit.aes = FALSE, vjust = 1,
            size = 3.6, colour = "grey30", lineheight = 0.9) +
  geom_step(linewidth = 0.8) +
  scale_colour_manual(values = cmap_set, guide = "none") +
  scale_x_log10(labels = label_comma(), breaks = 10^(0:7)) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 15)) +
  labs(x = "NTC cells needed for 10 linked molecules", y = "Percentage of sites (%)") +
  theme_manuscript(base_size = 13) + grid_both

FLOOR <- 1e-6
cmap_lev <- setNames(c("#08306B", "#2171B5", "#6BAED6", "#D94801"), LEVELS)
pd_d[, y := pmax(pct, FLOOR)]
pd <- ggplot(pd_d, aes(K, y, colour = Level, group = Level)) +
  geom_hline(yintercept = FLOOR, colour = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(data = pd_d[level == "cell"], aes(label = lab), vjust = -0.9, size = 3.4, show.legend = FALSE) +
  facet_wrap(~Set, nrow = 1) +
  scale_colour_manual(values = cmap_lev, name = NULL) +
  scale_y_log10(limits = c(FLOOR, 100), breaks = 10^(-6:2),
                labels = function(x) ifelse(x <= FLOOR, "0", sub("\\.?0+$", "", formatC(x, format = "f", digits = 6)))) +
  labs(x = "Linked molecules required on each allele of a heterozygous SNV",
       y = "Edited site \u00d7 unit pairs testable (%)") +
  theme_manuscript(base_size = 13, axis_text_size = 11) + grid_y +
  theme(legend.position = "bottom")

top <- patchwork::wrap_plots(pa, pb, pc, nrow = 1, widths = c(0.75, 1.6, 1.05)) +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
bottom <- patchwork::wrap_plots(pd, patchwork::plot_spacer(), nrow = 1, widths = c(2.2, 1))
fig <- patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(1, 0.95)) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 18))
save_figure_sd(fig, "TERATOMA_v8_linkage_snp_cells_allpanels", width = 19, height = 12.5,
               parts = c("a_sites", "b_pairs", "c_cells", "c_refs", "d_ased"), grid = "y", script = SCRIPT,
               src_stem = STEM)

ZERO <- 3e-5
STEPS <- levels(fn$Step)
fn[, `:=`(y = pmax(pct, ZERO), xi = as.integer(Level) + (step - 2) * 0.22,
          lab = ifelse(per_unit >= 10, format(round(per_unit), big.mark = ",", trim = TRUE),
                       format(signif(per_unit, 2), drop0trailing = TRUE, scientific = FALSE)))]
X_LAB <- sub(" ", "\n", LEVELS)
fn[, `:=`(hj = fifelse(level == "cell" & step == 2L, -0.3, 0.5), vj = fifelse(level == "cell" & step == 2L, 0.5, -1))]
pf <- ggplot(fn, aes(xi, y)) +
  geom_hline(yintercept = ZERO, colour = "grey75", linewidth = 0.4) +
  geom_line(aes(group = Level), colour = "grey60", linewidth = 0.7) +
  geom_point(aes(colour = Step), size = 3.4) +
  geom_text(aes(label = lab, colour = Step, hjust = hj, vjust = vj), size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c("#2171B5", "#D94801", "#7A0177"), STEPS), name = NULL,
                      guide = guide_legend(ncol = 1)) +
  scale_x_continuous(breaks = seq_along(LEVELS), labels = X_LAB, expand = expansion(add = 0.45)) +
  scale_y_log10(limits = c(ZERO, 150), breaks = c(ZERO, 10^(-4:1)),
                labels = function(x) ifelse(x <= ZERO, "0", sub("\\.?0+$", "", formatC(x, format = "f", digits = 4)))) +
  labs(x = NULL, y = "Percentage of truth site \u00d7 unit pairs (%)") +
  theme_manuscript(base_size = 13, axis_text_size = 11) + grid_y + theme(legend.position = "bottom")
save_figure_sd(pf, "TERATOMA_v8_linkage_snp_cells", width = 7, height = 6.5, parts = "ko_steps", grid = "y",
               script = SCRIPT, src_stem = STEM)
