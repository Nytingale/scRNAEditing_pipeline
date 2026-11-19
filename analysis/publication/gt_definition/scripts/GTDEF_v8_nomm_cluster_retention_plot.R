#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_nomm_cluster_retention_plot.R   (2026-10-18; draws from source data only since 2026-10-25)
# Response figure 8 (R2 major 3): does the discordant-mismatch filter (window 1) discard hyper-edited clusters?
#   a  A-to-G vs the other 11 mismatch types kept by each filter (window 2 = clustered-sites filter), alone and together
#   b  KO-responsive A-to-G sites kept by the discordant-mismatch filter, by the number of other A-to-G sites
#      within the truth w1 of the gate (nomm_w1_grid.tsv); facets gate x t
# Input : source_data/pm_gt/GTDEF_v8_nomm_cluster_retention.{a_bins,b_kept}.tsv
#         (GTDEF_v8_nomm_cluster_retention_source_data.R)
# Output: figures/pm_gt/GTDEF_v8_nomm_cluster_retention.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_plot_common.R"))
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/gt_definition/scripts/GTDEF_v8_nomm_cluster_retention_plot.R"
STEM <- "GTDEF_v8_nomm_cluster_retention"

b <- read_part(STEM, "a_bins")
m <- read_part(STEM, "b_kept")
levels(m$st) <- c("No local\nmismatch filters", "+ Discordant-\nmismatch filter", "+ Clustered-\nsites filter",
                  "+ Both filters")[match(levels(m$st), c("No local\ncontext filters", "+ Remove if near\nother MM types",
                                                          "+ Clustered", "+ Remove if near\nother MM types\n+ Clustered"))]
stopifnot(!anyNA(levels(m$st)))
levels(m$cls)[2] <- "Other 10 mismatch types with\ndecreased mismatch rates in KO"
BS <- 17
GT_V9 <- Sys.getenv("GT_V9") == "1"   # v9 truth labels (A-to-I = A-to-G + T-to-C)

pa <- ggplot(b, aes(clf, pct)) +
  geom_col(width = 0.5, fill = "#FDAE6B", colour = "grey20", linewidth = 0.2) +
  facet_grid(g ~ tl) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Number of sites with the same mismatch type kept",
       y = if (GT_V9) "KO-responsive A-to-G and\nT-to-C sites kept (%)" else "KO-responsive A-to-G sites kept (%)", tag = "a") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y +
  theme(panel.spacing.y = unit(2, "lines"))

cmap_cls <- setNames(c(cmap_varType[["AG"]], cmap_varType[["AC"]]), levels(m$cls))
pb <- ggplot(m, aes(st, pct, fill = cls)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "grey20", linewidth = 0.2) +
  facet_grid(g ~ tl) +
  scale_fill_manual(values = cmap_cls, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Sites kept (%)", tag = "b") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"), panel.spacing.y = unit(2, "lines"),
        legend.text = element_text(size = BS), legend.key.size = unit(1.4, "lines"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, lineheight = 0.85))

p <- (pb + labs(tag = "a")) / (pa + labs(tag = "b")) + plot_layout(heights = c(1.1, 1))
save_figure_sd(p, STEM, width = 12, height = 14.5, parts = c("a_bins", "b_kept"), script = SCRIPT)
ggsave(file.path(REPO_DIR, "review/response_figures_2026-10-15", paste0(STEM, ".png")), p,
       width = 12, height = 14.5, dpi = 150)
