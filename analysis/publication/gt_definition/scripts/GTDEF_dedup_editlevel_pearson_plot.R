#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_dedup_editlevel_pearson_plot.R   (2026-11-19)
# Response figure 15: Pearson's r of candidate-site A-to-G editing levels between two de-duplication methods,
# bars per tier (x) and method pair (fill), one panel per arm; r printed on each bar. Replaces the hex panels
# (GTDEF_hex_dedup_coord_vs_umi_plot.R) in the response. Fills: the three greens of the deduplication methods in
# cmap_proc (configuration tiles), in the order coordinate / umi_tools / UMI-aware MarkDuplicates.
# Input : source_data/pm_gt/GTDEF_dedup_editlevel_pearson.cor.tsv (GTDEF_dedup_editlevel_pearson_source_data.R)
# Output: figures/pm_gt/GTDEF_dedup_editlevel_pearson.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
SCRIPT <- "publication/gt_definition/scripts/GTDEF_dedup_editlevel_pearson_plot.R"
STEM <- "GTDEF_dedup_editlevel_pearson"

cr <- read_part(STEM, "cor")
cr[, pair := gsub("\\n", "\n", pair, fixed = TRUE)]
PAIR_LAB <- c("x: coordinate\ny: umi_tools" = "Coordinate-based vs umi_tools",
              "x: coordinate\ny: UMI-aware\nMarkDuplicates" = "Coordinate-based vs UMI-aware MarkDuplicates",
              "x: umi_tools\ny: UMI-aware\nMarkDuplicates" = "umi_tools vs UMI-aware MarkDuplicates")
stopifnot(all(cr$pair %in% names(PAIR_LAB)))
cr[, Pair := factor(PAIR_LAB[pair], PAIR_LAB)]
cr[, tier := factor(tier, c("Arm", "Sample", "Cell type"))]
cr[, arm := factor(arm, c("NTC", "ADAR1 KO", "ADARB1 KO"))]
cr[, `:=`(lab = ifelse(is.na(r), "not run", sprintf("%.3f", r)), r0 = fifelse(is.na(r), 0, r),
          ty = fifelse(is.na(r), 0.02, r - 0.02), hj = fifelse(is.na(r), 0, 1),
          txt_col = fifelse(is.na(r), "grey30", fifelse(Pair == PAIR_LAB[[1]], "white", "black")))]
BS <- 13

p <- ggplot(cr, aes(tier, r0, fill = Pair)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(y = ty, label = lab, colour = txt_col, hjust = hj, group = Pair),
            position = position_dodge(width = 0.85), angle = 90, size = 3.8, show.legend = FALSE) +
  scale_colour_identity() +
  facet_grid(. ~ arm) +
  scale_fill_manual(values = setNames(cmap_proc[c("De-duplicated (coordinate)", "De-duplicated (umi_tools)",
                                                  "De-duplicated (UMI-aware MarkDuplicates)")], PAIR_LAB), name = NULL) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Pearson's r of A-to-G editing levels") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y +
  theme(legend.position = "bottom", legend.text = element_text(size = BS), panel.spacing.x = unit(10, "pt"))
save_figure_sd(p, STEM, width = 8, height = 5, parts = "cor", script = SCRIPT)
