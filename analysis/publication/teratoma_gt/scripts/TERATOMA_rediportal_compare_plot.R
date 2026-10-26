#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_rediportal_compare_plot.R   (HANDOFF 2026-09-29 UPDATE-25)
# Ranking of the reported configurations (full benchmark call set) under the v6 KO-depletion
# truth vs under REDIportal v2 as truth (gt_score_rediportal_v6.R). Same rank-sum as the
# cross-tier grid: F1 rank within tier, summed over the three tiers.
#
# Inputs
#   human/pm/gt_eval_v6/scoring/sweep_infneg/teratoma_inputs/gt_tier_scoring.tsv   (v6, t = 50%)
#   human/pm/gt_eval_v6/rediportal/scoring/gt_tier_scoring_{cov10_alt1,alt3,alt3_e05}.tsv (2026-10-01)
# Outputs (figures/pm_gt/)
#   TERATOMA_rediportal_rank_vs_v6.svg      overall rank, ground truth vs each REDIportal variant
#   TERATOMA_rediportal_f1_by_tier.svg      per-tier F1, ground truth vs REDIportal (>= 3 G, edit >= 5%)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_rediportal_compare_plot.R"

RT <- NINE_TRUTHS[startsWith(truth, "redi_")]
VAR_LAB <- setNames(RT$label, sub("^redi_", "", RT$truth))
F1_VAR <- "alt3_e05"
TIER_LAB <- c(dataset = "Arm", sample = "Sample", celltype = "Cell type")
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
SRC <- c(v6 = file.path(SCORING_DIR, "gt_tier_scoring.tsv"), setNames(RT$file, names(VAR_LAB)))

rank_sum <- function(s) {
  s <- copy(s)[, f1 := fifelse(precision + recall_vs_union > 0,
                               2 * precision * recall_vs_union / (precision + recall_vs_union), NA_real_)]
  s[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
  r <- s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
  r[, overall_rank := frank(list(rank_sum, -mean_f1), ties.method = "min")]   # ties: higher mean F1
  list(tier = s, overall = r)
}
sc <- lapply(SRC, function(f) rank_sum(fread(f)))

# ---- (a) overall rank: v6 vs REDIportal ------------------------------------------
ov <- rbindlist(lapply(names(VAR_LAB), function(v)
  merge(sc$v6$overall[, c(COMBO, "overall_rank"), with = FALSE],
        sc[[v]]$overall[, c(COMBO, "overall_rank"), with = FALSE],
        by = COMBO, suffixes = c("_v6", "_redi"))[, variant := v]))
ov[, `:=`(Variant = factor(VAR_LAB[variant], VAR_LAB), method = std_caller(caller))]
rho <- ov[, .(rho = cor(overall_rank_v6, overall_rank_redi, method = "spearman"),
              top10_shared = sum(overall_rank_v6 <= 10 & overall_rank_redi <= 10)), by = Variant]
rho[, lab := sprintf("Spearman's rho = %.2f\ntop-10 shared: %d", rho, top10_shared)]
fwrite(ov, file.path(SRC_DIR, "TERATOMA_rediportal_rank_vs_v6.tsv"), sep = "\t")

STEM <- "TERATOMA_rediportal_rank_vs_v6"
p <- ggplot(ov, aes(overall_rank_v6, overall_rank_redi)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = method, shape = strand_mode), size = 2, alpha = 0.85) +
  geom_text(data = rho, aes(x = 1, y = Inf, label = lab), hjust = 0, vjust = 1.2, size = 3.2) +
  facet_wrap(~ Variant, nrow = 1) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) + scale_y_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
  scale_colour_manual_cvd(values = cmap_method, name = NULL) +
  scale_shape_manual(values = c(annotation = 16, hybrid = 17, caller = 15), labels = label_strandMode, name = NULL) +
  labs(x = "Overall rank, ground truth (>= 50% depletion)", y = "Overall rank, REDIportal as truth",
       title = "Cross-tier rank of each configuration: ground truth vs REDIportal") +
  theme_manuscript(facets = TRUE) + grid_both +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.box = "vertical")
save_figure(p, STEM, width = 12, height = 5.8, sources = SRC, data = ov, grid = "both", script = SCRIPT)

# ---- (b) per-tier F1: ground truth vs REDIportal F1_VAR (square panels, equal x/y range per tier) --------------------------------------
tf <- merge(sc$v6$tier[, c("tier", COMBO, "f1"), with = FALSE],
            sc[[F1_VAR]]$tier[, c("tier", COMBO, "f1", "precision", "recall_vs_union"), with = FALSE],
            by = c("tier", COMBO), suffixes = c("_v6", "_redi"))
tf[, `:=`(Tier = factor(TIER_LAB[tier], TIER_LAB), method = std_caller(caller))]
rt <- tf[, .(lab = sprintf("Spearman's rho = %.2f", cor(f1_v6, f1_redi, method = "spearman"))), by = Tier]
lims <- tf[, .(lo = min(f1_v6, f1_redi, na.rm = TRUE), hi = max(f1_v6, f1_redi, na.rm = TRUE)), by = Tier]
lims <- lims[, .(f1_v6 = c(lo, hi), f1_redi = c(lo, hi)), by = Tier]
fwrite(tf, file.path(SRC_DIR, "TERATOMA_rediportal_f1_by_tier.tsv"), sep = "\t")
STEM <- "TERATOMA_rediportal_f1_by_tier"
p <- ggplot(tf, aes(f1_v6, f1_redi)) +
  geom_blank(data = lims) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = method, shape = strand_mode), size = 2, alpha = 0.85) +
  geom_text(data = rt, aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.3, size = 3.2) +
  facet_wrap(~ Tier, scales = "free", nrow = 1) +
  scale_colour_manual_cvd(values = cmap_method, name = NULL) +
  scale_shape_manual(values = c(annotation = 16, hybrid = 17, caller = 15), labels = label_strandMode, name = NULL) +
  labs(x = "F1, ground truth (>= 50% depletion)", y = paste0("F1, ", VAR_LAB[[F1_VAR]]),
       title = "Per-tier F1 under the two truths") +
  theme_manuscript(facets = TRUE) + grid_both +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.box = "vertical")
save_figure(p, STEM, width = 11, height = 5.5, sources = SRC[c("v6", F1_VAR)], data = tf, grid = "both", script = SCRIPT)

# ---- summary ------------------------------------------------------------------------
print(rho[, .(Variant, rho = round(rho, 3), top10_shared)])
print(rt)
top <- function(x, k = 10) x$overall[order(overall_rank)][1:k, c("overall_rank", COMBO, "rank_sum"), with = FALSE]
cat("\nv6 top 10\n"); print(top(sc$v6))
for (v in names(VAR_LAB)) { cat("\n", v, " top 10\n", sep = ""); print(top(sc[[v]])) }
best_caller <- rbindlist(lapply(names(SRC), function(s)
  sc[[s]]$overall[, .(best_rank = min(overall_rank)), by = caller][, truth := s]))
cat("\nbest overall rank per caller\n"); print(dcast(best_caller, caller ~ truth, value.var = "best_rank"))
