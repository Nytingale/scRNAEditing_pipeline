#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_split_cells_plot.R   (2026-10-14)
# Supplementary split-by-cells test: NTC cells split at random 50/50 within teratoma x cell type (seed 20261013).
# Calls from half B scored (all three tiers) against (top) the full v8 truth and (bottom) a truth rebuilt from
# half A with every rule frozen (reads disjoint from the calls). x = overall cross-tier rank of the 58 reported
# configurations in the main analysis (all cells), y = rank on half B. Top row = depth effect; bottom vs top =
# shared-reads effect at equal depth.
#
# Inputs : human/pm/gt_eval_v8/split_cells/compare_config_rank.tsv, truth_overlap.tsv
#          (human/scripts/gt_split_compare.R, human/pm/scripts/run_gt_split_score.sh)
# Output : figures/pm_gt/TERATOMA_v8_scatter_rank_split_cells.svg
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_split_cells_plot.R"

SDIR <- file.path(V8_DIR, "split_cells")
SRC <- file.path(SDIR, c("compare_config_rank.tsv", "truth_overlap.tsv"))
stopifnot(all(file.exists(SRC)))
BS <- 14

x <- fread(SRC[1])
stopifnot(x[, .N, by = .(gate, t)][, all(N == 58L)])
d <- melt(x, id.vars = c("gate", "t", "config", "caller", "bam_type", "rank_main"),
          measure.vars = c("rank_B_full", "rank_B_A"), variable.name = "truth", value.name = "rank_split")
d[, truth := sub("^rank_", "", truth)]
ov <- fread(SRC[2])

st <- d[, .(rho = cor(rank_main, rank_split, method = "spearman"),
            top1_split = config[rank_split == 1][1]), by = .(gate, t, truth)]
st <- merge(st, x[, .(rho_AB = cor(rank_B_full, rank_B_A, method = "spearman")), by = .(gate, t)], by = c("gate", "t"))
st <- merge(st, ov[, .(gate, t, pos_A, pos_full)], by = c("gate", "t"))
print(st[order(truth, gate, t)])
fmt <- function(n) format(n, big.mark = ",", trim = TRUE)
st[, lab := ifelse(truth == "B_full",
                   sprintf("\u03C1 vs all cells = %.3f\n%s positives", rho, fmt(pos_full)),
                   sprintf("\u03C1 vs all cells = %.3f\n\u03C1 vs top row = %.3f\n%s positives", rho, rho_AB, fmt(pos_A)))]

T_LAB <- c(t25 = "\u2265 25% editing reduction in KO", t50 = "\u2265 50% editing reduction in KO", t75 = "\u2265 75% editing reduction in KO")
ROW_LAB <- c(B_full = "Half-B calls,\ntruth from all cells", B_A = "Half-B calls,\ntruth from half A")
fac <- function(z) {
  z[, `:=`(col = factor(sprintf("%s\n%s", V8_GATE_LAB2[gate], T_LAB[t]),
                        levels = as.vector(outer(V8_GATE_LAB2, T_LAB, function(g, tt) sprintf("%s\n%s", g, tt)))[c(1, 3, 5, 2, 4, 6)]),
           row = factor(ROW_LAB[truth], levels = ROW_LAB))]
}
fac(d); fac(st)
d[, `:=`(method = caller_factor(caller), proc = relabel(bam_type, label_bamType))]
PROC_LEVELS <- intersect(label_bamType[c("sampleRmDup", "sampleUmiDup", "recal", "recalq20", "recaldefault")], d$proc)
d[, proc := factor(proc, PROC_LEVELS)]
shape_proc <- setNames(c(21, 24, 22, 23, 25)[seq_along(PROC_LEVELS)], PROC_LEVELS)
assert_cvd_safe(cmap_method, "scale_fill_manual")

# both axes reversed: rank 1 at the top right
p <- ggplot(d, aes(rank_main, rank_split)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(fill = method, shape = proc), colour = "grey15", stroke = 0.3, size = 3.4, alpha = 0.7) +
  geom_text(data = st, aes(x = 58, y = 1, label = lab), hjust = 0, vjust = 1, size = BS / 3.3, lineheight = 0.95,
            inherit.aes = FALSE) +
  facet_grid(row ~ col) +
  scale_x_reverse(breaks = c(58, 40, 20, 1), expand = expansion(add = 1.5)) +
  scale_y_reverse(breaks = c(1, 20, 40, 58), expand = expansion(add = 1.5)) +
  scale_fill_manual(values = cmap_method, name = NULL,
                    guide = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, size = 4.5, alpha = 0.9))) +
  scale_shape_manual(values = shape_proc, name = NULL,
                     guide = guide_legend(nrow = 1, order = 2, override.aes = list(fill = "grey60", size = 4.5, alpha = 0.9))) +
  labs(x = "Overall cross-tier rank, all NTC cells", y = "Overall cross-tier rank, half B of the NTC cells",
       title = "Ranking with calls and ground truth from disjoint halves of the NTC cells",
       caption = "\u03C1: Spearman's \u03C1 of the ranks") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_both +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        legend.box = "vertical", legend.spacing.y = unit(2, "pt"))

save_figure(p, "TERATOMA_v8_scatter_rank_split_cells", width = 20, height = 9.5, sources = SRC,
            data = d[, .(gate, t, truth, config, caller, bam_type, rank_main, rank_split)], grid = "both", script = SCRIPT)
