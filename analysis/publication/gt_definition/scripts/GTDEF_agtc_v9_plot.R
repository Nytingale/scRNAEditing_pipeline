#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_agtc_v9_plot.R
#
# Draws the source data of GTDEF_agtc_v9_source_data.R:
#   a  per mismatch type (read orientation), % of gate sites with editing reduction in KO >= t;
#      rows = dataset, columns = gate
#   b  where the gene bodies lie relative to the read strand, by KO-reduction bin, for A-to-G / T-to-C / other
#
# Input : source_data/pm_gt/GTDEF_agtc_v9.{pass,gene}.tsv
# Output: figures/pm_gt/GTDEF_panel_agtc_v9_response.{svg,pdf}
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/gt_definition/scripts/GTDEF_agtc_v9_plot.R"
STEM <- "GTDEF_agtc_v9"
BS <- 12
DS <- c("Teratoma", "Aorta")
GATES <- c("Power-gated", "Depth-gated")
GRP_LAB <- c(AG = "A-to-G", TC = "T-to-C", other = "Other 10 types")
cmap_grp <- c("A-to-G" = "#D55E00", "T-to-C" = "#0072B2", "Other 10 types" = "grey60")
GENE_LEV <- c("Gene on the opposite strand", "Gene on the read strand", "Genes on both strands", "No gene")
cmap_gene <- setNames(c("#009E73", "#CC79A7", "#56B4E9", "#E69F00"), GENE_LEV)
BIN_LAB <- c("< 0", "0-25", "25-50", "50-75", "\u2265 75")

pass <- read_part(STEM, "pass")
pass[, `:=`(dataset = factor(dataset, DS), gate = factor(gate, GATES),
            grp = factor(GRP_LAB[fifelse(type %in% c("AG", "TC"), type, "other")], GRP_LAB))]
pa <- ggplot(pass, aes(100 * t, pct, group = type, colour = grp)) +
  geom_line(data = pass[grp == GRP_LAB[["other"]]], linewidth = 0.5, alpha = 0.8) +
  geom_point(data = pass[grp == GRP_LAB[["other"]]], size = 1.1, alpha = 0.8) +
  geom_line(data = pass[grp != GRP_LAB[["other"]]], linewidth = 1.2) +
  geom_point(data = pass[grp != GRP_LAB[["other"]]], size = 2) +
  facet_grid(dataset ~ gate) +
  scale_colour_manual(values = cmap_grp, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Editing reduction in KO threshold, t (%)", y = "Sites with \u2265 t% editing\nreduction in KO (%)") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_both +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

gene <- read_part(STEM, "gene")
gene[, `:=`(dataset = factor(dataset, DS), gate = factor(gate, GATES), grp = factor(GRP_LAB[grp], GRP_LAB),
            bin = factor(bin, BIN_LAB), gene_class = factor(gene_class, GENE_LEV))]
pb <- ggplot(gene[n_bin > 0], aes(bin, pct, fill = gene_class)) +
  geom_col(width = 0.85, colour = "grey20", linewidth = 0.15) +
  facet_grid(dataset ~ grp + gate) +
  scale_fill_manual(values = cmap_gene, name = NULL, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 100, 25)) +
  labs(x = "Editing reduction in KO (%)", y = "Sites (%)") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y +
  theme(legend.position = "bottom", panel.spacing = unit(0.6, "lines"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

p <- (pa / pb) + plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 16, face = "bold"))
save_figure_sd(p, "GTDEF_panel_agtc_v9_response", 11, 10.5, parts = c("pass", "gene"), grid = "y",
               script = SCRIPT, src_stem = STEM)
ggsave(file.path(FIG_DIR, "GTDEF_panel_agtc_v9_response.pdf"), p, width = 11, height = 10.5, device = cairo_pdf)
message("Done")
