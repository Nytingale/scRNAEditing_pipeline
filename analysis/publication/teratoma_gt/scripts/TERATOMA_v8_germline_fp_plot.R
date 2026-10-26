#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_germline_fp_plot.R   (2026-10-16)
# Supplementary, truth-independent check: share of each configuration's discovered A-to-G sites that are PGP1
# germline A>G SNVs (gt_score_germline_v8.R germ_frac: per unit germline calls / A-to-G calls, mean over units,
# then over the arm, sample and cell-type tiers), stacked by germline genotype, for the 58 ranked configurations.
# Rows: configurations (tiles as TERATOMA_v8_rank_robustness_all_tiles), ordered by the overall cross-tier rank
# (rank-sum rank over the nine v8 truths, as TERATOMA_v8_window_state_rank).
#
# Inputs : human/pm/gt_eval_v8/state_summary/config_ranks.tsv (human/scripts/gt_v8_state_caller_summary.R)
#          human/pm/gt_eval_v8/scoring/mismatch/config_rank.tsv (per KO truth ranks, for the rho table)
#          human/pm/gt_eval_v8/germline/scoring/mismatch/<tier>/germline_scores.tsv, germline/germline_counts.tsv
#          (human/scripts/gt_v8_germline_sites.R, human/pm/scripts/run_gt_v8_germline.sh)
# Output : figures/pm_gt/TERATOMA_v8_scatter_germline_fp_vs_rank.svg (stem kept for the response build)
#          human/pm/gt_eval_v8/germline/germline_vs_rank.tsv, germline_best_per_caller.tsv
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
suppressPackageStartupMessages(library(patchwork))
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_germline_fp_plot.R"

GDIR <- file.path(V8_DIR, "germline")
TIERS <- c("arm", "sample", "celltype")
GFILES <- file.path(GDIR, "scoring", V8_STATE, TIERS, "germline_scores.tsv")
RANKS <- file.path(V8_DIR, "state_summary", "config_ranks.tsv")
SRC <- c(RANKS, v8_config_rank(), GFILES, file.path(GDIR, "germline_counts.tsv"))
stopifnot(all(file.exists(SRC)))
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
BASE <- 16

r <- fread(RANKS)[state == V8_STATE]
stopifnot(setequal(r$truth, V8_TRUTH_IDS), r[, .N, by = truth][, all(N == 58L)])
rs <- r[, .(rank_sum = sum(overall_rank), mean_f1 = mean(mean_f1)), by = c("config", COMBO)]
rs[, rs_rank := frank(list(-mean_f1, rank_sum), ties.method = "min")]

GT <- c(germ_frac_hom = "1/1", germ_frac_gt01 = "0/1", germ_frac_gtdot1 = "./1", germ_frac_gt1 = "1")
MET <- c("germ_frac", names(GT))
g <- rbindlist(lapply(GFILES, fread))[config %in% rs$config]
stopifnot(g[, .N, by = tier][, all(N == 58L)], all(MET %in% names(g)))
gw <- dcast(g, config ~ tier, value.var = MET)
for (m in MET) gw[, (m) := rowMeans(.SD), .SDcols = paste(m, TIERS, sep = "_")]
gw <- gw[, c("config", MET), with = FALSE]
stopifnot(gw[, all(abs(germ_frac - Reduce(`+`, .SD)) < 1e-9), .SDcols = names(GT)])
d <- merge(rs, gw, by = "config")
stopifnot(nrow(d) == 58L)
setorder(d, rs_rank, -mean_f1)
fwrite(d, file.path(GDIR, "germline_vs_rank.tsv"), sep = "\t")

best <- d[order(rs_rank), .SD[1], by = caller][, .(caller, config, rs_rank, germ_frac, germ_frac_hom)]
fwrite(best, file.path(GDIR, "germline_best_per_caller.tsv"), sep = "\t")
cat("Best-ranked configuration per caller (germline share, %):\n")
print(best[, .(caller, config, rs_rank, pct = round(100 * germ_frac, 2), hom_pct = round(100 * germ_frac_hom, 2))])
cat("\nLowest germline share (%):\n"); print(d[order(germ_frac)][1:5, .(config, rs_rank, pct = round(100 * germ_frac, 3))])
cat("\nHighest germline share (%):\n"); print(d[order(-germ_frac)][1:5, .(config, rs_rank, pct = round(100 * germ_frac, 3))])
cat("\nSpearman's rho, overall rank vs germline share:", round(cor(d$rs_rank, d$germ_frac, method = "spearman"), 3), "\n")
ko <- fread(v8_config_rank())[, .(gate, t, config, rank_main = overall_rank)]
kr <- merge(ko, d[, c("config", MET), with = FALSE], by = "config")
print(kr[, .(rho_all = round(cor(rank_main, germ_frac, method = "spearman"), 2),
             rho_hom = round(cor(rank_main, germ_frac_hom, method = "spearman"), 2)), by = .(gate, t)])

cnt <- fread(file.path(GDIR, "germline_counts.tsv"))
n_gt <- setNames(vapply(GT, function(x) cnt[step == paste0("kept, GT ", x), n], integer(1)), GT)
stopifnot(sum(n_gt) == cnt[step == "not multiallelic", n])
GT_LAB <- setNames(sprintf(c("Homozygous, GT 1/1 (n = %d)", "Heterozygous, GT 0/1 (n = %d)",
                             "Heterozygous, GT ./1 (n = %d)", "Hemizygous chrX, GT 1 (n = %d)"), n_gt), names(GT))
cmap_gt <- setNames(c("#08519C", "#6BAED6", "#C6DBEF", "#E6550D"), GT_LAB)

N <- nrow(d)
d[, `:=`(y = N + 1L - .I, Preprocessing = relabel(bam_type, label_bamType), Caller = relabel(caller, label_caller),
         Stranding = relabel(strand_mode, label_strandMode), `SNP removal` = relabel(snp_subtract, label_snpSubtract))]
long <- melt(d, id.vars = c("config", "y", "rs_rank"), measure.vars = names(GT), variable.name = "gt", value.name = "frac")
long[, Genotype := factor(GT_LAB[as.character(gt)], rev(GT_LAB))]

CMAPS <- list(Preprocessing = cmap_proc, Caller = cmap_method, Stranding = cmap_strand, `SNP removal` = cmap_snp)
Y_LIM <- c(0.5, N + 0.5)
tile_col <- function(col) {
  td <- d[, .(y, value = get(col), column = col)]
  ggplot(td, aes(column, y, fill = value)) +
    geom_tile(width = 1, height = 0.78) +
    scale_fill_manual(values = CMAPS[[col]], name = col,
                      breaks = if (col == "Caller") levels(caller_factor(d$caller))[
                                 levels(caller_factor(d$caller)) %in% td$value]
                               else intersect(names(CMAPS[[col]]), unique(td$value))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 1, title.position = "top")) +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = BASE) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
          axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(5, 0, 5, 0))
}
p_rank <- ggplot(d, aes("Rank", y, label = rs_rank)) +
  geom_text(size = BASE / 3.2) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_manuscript(base_size = BASE) +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        plot.margin = margin(5, 2, 5, 0))
p_bar <- ggplot(long, aes(100 * frac, y, fill = Genotype)) +
  geom_col(orientation = "y", width = 0.78, colour = NA) +
  scale_fill_manual(values = cmap_gt, breaks = unname(GT_LAB), name = "Germline genotype",
                    guide = guide_legend(ncol = 1, title.position = "top")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
  labs(x = "Germline SNVs among discovered\nA-to-G sites (%, mean over tiers)", y = NULL) +
  theme_manuscript(base_size = BASE) + grid_x +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey77"), plot.margin = margin(5, 10, 5, 4))

TILE_IN <- 0.36
panels <- c(list(p_rank), lapply(names(CMAPS), tile_col), list(p_bar))
p <- wrap_plots(panels, nrow = 1,
                widths = unit(c(0.45, rep(TILE_IN, length(CMAPS)), 1), c("in", rep("in", length(CMAPS)), "null"))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.justification = "center", legend.title = element_text(face = "bold"))
save_figure(p, "TERATOMA_v8_scatter_germline_fp_vs_rank", width = 14, height = 2.6 + 0.27 * N, sources = SRC,
            data = long[, .(config, rs_rank, gt, frac)], grid = "x", script = SCRIPT)
