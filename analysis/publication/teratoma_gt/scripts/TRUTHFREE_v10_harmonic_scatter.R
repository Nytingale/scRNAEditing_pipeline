#!/usr/bin/env Rscript
# TRUTHFREE_v10_harmonic_scatter.R   (user 2026-11-14: Response Figure 12 on the 76 configurations)
# tmp/meeting_v9_2026-11-09/truthfree_scatter_harmonic.R (RESPONSE=1) for the 76 ranked configurations (58 published,
# 2 Cellsnp-lite hybrid, 16 UMI-aware MarkDuplicates; umi_tools seeded): truth-free measures (A-to-G percentage, number
# of A-to-G sites, harmonic score) against mean F1, precision and recall, teratoma (a) and held-out aorta (b), with
# least-squares lines and Spearman's rho. Also prints the statistics quoted in the response (Reviewer 1, major 1):
# site-count vs percentage ranks in the teratoma arm, ranges and rho of the KO measures, harmonic vs average-rank rho.
# Harmonic score = 2ab / (a + b), a = A-to-G sites / max, b = A-to-G percentage / max over the 76 configurations.
# Input : human/pm/gt_eval_v9/configs76/configs76_seeded.tsv (human/scripts/gt_v10_configs76.R)
# Output: review/response_figures_2026-10-15/TRUTHFREE_v10_harmonic_teratoma_aorta.png,
#         publication/teratoma_gt/figures/pm_gt76/TRUTHFREE_v10_harmonic_teratoma_aorta.pdf,
#         publication/teratoma_gt/source_data/pm_gt76/TRUTHFREE_v10_harmonic_scatter.tsv
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
setwd(Sys.getenv("SCRNA_DIR", "."))
source("plotting/plotting_functions.R")
OUT <- "review/response_figures_2026-10-15"
SD <- "publication/teratoma_gt/source_data/pm_gt76"
FIG <- "publication/teratoma_gt/figures/pm_gt76"
for (o in c(SD, FIG)) dir.create(o, recursive = TRUE, showWarnings = FALSE)
BASE <- 10
N_CONFIG <- 76L
CALLERS <- c(bcftools = "BCFtools", cellsnp = "Cellsnp-lite", gatk = "GATK", jacusa2 = "JACUSA2",
             reditools3 = "REDItools3", scallele = "scAllele")
col_caller <- setNames(cmap_method[names(CALLERS)], CALLERS)
MEAS <- c(pct = "A-to-G percentage", sites = "Sites (thousands)", harm = "Harmonic score")
MET <- c(f1 = "Mean F1", precision = "Precision", recall = "Recall")
BAMS <- c(sampleRmDup = "De-duplicated (coordinate)", sampleUmiDup = "De-duplicated (umi_tools)",
          sampleUmiMd = "De-duplicated (UMI-aware MarkDuplicates)", recal = "Recalibrated",
          recalq20 = "Recalibrated (Q20)", recaldefault = "Recalibrated (default)")
shp_proc <- setNames(c(21, 22, 23, 24, 25, 4), BAMS)

x0 <- fread("human/pm/gt_eval_v9/configs76/configs76_seeded.tsv")
stopifnot(nrow(x0) == N_CONFIG, !anyNA(x0))
x0[, bam := tstrsplit(config, "__")[[2]]]
stopifnot(all(x0$bam %in% names(BAMS)))
x <- rbind(
  x0[, .(dataset = "Teratoma (human KO)", config, caller, bam, sites = ag_te, pct = pct_te,
         f1 = f1_te, precision = precision_te, recall = recall_te)],
  x0[, .(dataset = "Aorta (mouse conditional KO)", config, caller, bam, sites = ag_ao, pct = pct_ao,
         f1 = f1_ao, precision = precision_ao, recall = recall_ao)])
x[, `:=`(Caller = factor(CALLERS[caller], CALLERS), Processing = factor(BAMS[bam], BAMS))]
x[, harm := { a <- sites / max(sites); b <- pct / max(pct); 2 * a * b / (a + b) }, by = dataset]
x[, `:=`(rank_sites = frank(-sites, ties.method = "min"), rank_pct = frank(-pct, ties.method = "min"),
         rank_f1 = frank(-f1, ties.method = "min")), by = dataset]
x[, avg_rank := (rank_sites + rank_pct) / 2]
fwrite(x, file.path(SD, "TRUTHFREE_v10_harmonic_scatter.tsv"), sep = "\t")

sp <- function(a, b) cor(a, b, method = "spearman")
cat("\n== Statistics for the response (", N_CONFIG, " configurations) ==\n", sep = "")
for (ds in unique(x$dataset)) {
  y <- x[dataset == ds]
  cat(sprintf("\n%s\n  precision %.2f-%.2f, recall %.2f-%.2f\n", ds, min(y$precision), max(y$precision),
              min(y$recall), max(y$recall)))
  cat(sprintf("  rho F1~recall %.2f, F1~precision %.2f\n", sp(y$f1, y$recall), sp(y$f1, y$precision)))
  cat(sprintf("  rho F1~sites %.2f, F1~pct %.2f, F1~harmonic %.2f, F1~(-average rank) %.2f\n",
              sp(y$f1, y$sites), sp(y$f1, y$pct), sp(y$f1, y$harm), sp(y$f1, -y$avg_rank)))
  cat(sprintf("  rho sites~pct %.2f\n", sp(y$sites, y$pct)))
  top <- y[which.max(harm)]
  cat(sprintf("  top harmonic: %s (F1 %.3f, rank %d by F1); best F1 %.3f (%s); gap %.4f\n", top$config, top$f1,
              top$rank_f1, max(y$f1), y[which.max(f1), config], max(y$f1) - top$f1))
}
te <- x[dataset == "Teratoma (human KO)"]
te[, move := rank_pct - rank_sites]
cat(sprintf("\nTeratoma arm: configurations moving >= 20 places between site-count and percentage ranks: %d of %d\n",
            te[abs(move) >= 20, .N], nrow(te)))
print(te[abs(move) >= 20][order(caller, move), .(config, rank_sites, rank_pct, move)])
print(te[abs(move) >= 20, .(n = .N, sites_ranks = paste(range(rank_sites), collapse = "-"),
                            pct_ranks = paste(range(rank_pct), collapse = "-")), by = .(caller, direction = sign(move))])

d <- melt(x, id.vars = c("dataset", "config", "Caller", "Processing"),
          measure.vars = list(xv = c("pct", "sites", "harm")), variable.name = "measure", value.name = "xv")
d[, measure := names(MEAS)[measure]]
d <- merge(d, x[, .(dataset, config, f1, precision, recall)], by = c("dataset", "config"))
d <- melt(d, measure.vars = names(MET), variable.name = "metric", value.name = "yv")
d[measure == "sites", xv := xv / 1e3]
d[, `:=`(Measure = factor(MEAS[measure], MEAS), Metric = factor(MET[as.character(metric)], MET))]
rho <- d[, .(rho = sp(xv, yv), p = spearman_p(xv, yv), n = .N), by = .(dataset, Measure, Metric)]
stopifnot(rho[, all(n == N_CONFIG)])
rho[, lab := sprintf("Spearman's \u03c1 = %.2f\n%s", rho, p_lab(p))]
print(dcast(rho, dataset + Metric ~ Measure, value.var = "rho"), digits = 2)

TITLE <- c("Teratoma (human KO)" = "Human teratoma", "Aorta (mouse conditional KO)" = "Mouse aorta")
plot_one <- function(ds) {
  dd <- d[dataset == ds]; r <- rho[dataset == ds]
  ggplot(dd, aes(xv, yv)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey30", linetype = "dashed", linewidth = 0.5) +
    geom_point(aes(fill = Caller, shape = Processing), size = 1.6, stroke = 0.3, colour = "grey15", alpha = 0.9) +
    geom_text(data = r, aes(x = -Inf, y = Inf, label = lab), hjust = -0.06, vjust = 1.2, size = 2.9, lineheight = 0.95, inherit.aes = FALSE) +
    facet_grid(Metric ~ Measure, scales = "free", switch = "both") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    scale_fill_manual(values = col_caller, name = NULL,
                      guide = guide_legend(override.aes = list(shape = 21, size = 3), nrow = 3, order = 1)) +
    scale_shape_manual(values = shp_proc, name = NULL, drop = FALSE,
                       guide = guide_legend(override.aes = list(fill = "grey60", size = 3), ncol = 1, order = 2)) +
    theme_manuscript(base_size = BASE, facets = TRUE) + grid_both +
    theme(strip.placement = "outside", legend.position = "bottom", legend.box = "vertical",
          legend.text = element_text(size = BASE), aspect.ratio = 1,
          strip.text = element_text(size = BASE), strip.background = element_blank(), panel.spacing.x = unit(26, "pt"),
          axis.text = element_text(size = BASE - 1), plot.title = element_text(size = BASE + 1),
          legend.box.spacing = unit(2, "pt")) +
    labs(x = NULL, y = NULL, title = TITLE[[ds]])
}
p <- (plot_one("Teratoma (human KO)") / plot_one("Aorta (mouse conditional KO)")) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "horizontal", plot.tag = element_text(size = BASE + 4, face = "bold"))
ggsave(file.path(OUT, "TRUTHFREE_v10_harmonic_teratoma_aorta.png"), p, width = 7.5, height = 14, dpi = 300, bg = "white")
ggsave(file.path(FIG, "TRUTHFREE_v10_harmonic_teratoma_aorta.pdf"), p, width = 7.5, height = 14, device = cairo_pdf)
