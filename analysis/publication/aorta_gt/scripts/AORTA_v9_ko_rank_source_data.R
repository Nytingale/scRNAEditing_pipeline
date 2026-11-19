#!/usr/bin/env Rscript
# ============================================================================
# AORTA_v9_ko_rank_source_data.R   (2026-11-19; response R1 major 3: held-out aorta ranking and cross-species agreement)
# The 76 configurations scored on the six v9 KO truth sets of the mouse aorta (arm and sample tiers) and the
# human teratoma (arm, sample and cell-type tiers):
#   ranks    aorta rank of every configuration under each KO truth set
#   configs  aorta configurations in figure order (worst rank, then mean rank), as response figure 7a
#   margin   per aorta truth set, BCFtools' best mean F1 minus the best mean F1 of any other caller
#   cross    mean F1 across the six KO truth sets on each dataset and the rank by it
# Input : mouse/aorta_gse255661/gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv
#         human/pm/gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv
# Output: source_data/aorta_gt/AORTA_v9_ko_rank.{ranks,configs,margin,cross}.tsv (AORTA_v9_ko_rank_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "AORTA_plot_common.R"))
STEM <- "AORTA_v9_ko_rank"
N_CFG <- 76L

label_bamType[c("sampleUmiDup", "sampleUmiMd")] <- c("De-duplicated (umi_tools)", "De-duplicated (UMI-aware MarkDuplicates)")
read_rank <- function(f) {
  r <- fread(f)[, truth := paste(gate, t, sep = "_")]
  stopifnot(setequal(r$truth, V8_KO_TRUTHS), r[, .N, by = truth][, all(N == N_CFG)])
  r[, c(COMBO, "config", "truth", "mean_f1", "rank_sum", "overall_rank"), with = FALSE]
}
ao <- read_rank(file.path(AO_DIR, "gt_eval_v9", "scoring_pwcov76", "mismatch", "config_rank.tsv"))
te <- read_rank(file.path(REPO_DIR, "human", "pm", "gt_eval_v9", "scoring_pwcov76", "mismatch", "config_rank.tsv"))
add_labels <- function(d) d[, `:=`(label = cfg_label(d), Preprocessing = relabel(bam_type, label_bamType),
                                   Caller = relabel(caller, label_caller), Stranding = relabel(strand_mode, label_strandMode),
                                   `SNP removal` = relabel(snp_subtract, label_snpSubtract))]

# ---- aorta ranks and figure order ----------------------------------------------------
cfg <- ao[, .(worst = max(overall_rank), best = min(overall_rank), mean_rank = mean(overall_rank)), by = COMBO]
setorder(cfg, worst, mean_rank)
cfg[, y := .N + 1L - .I]
add_labels(cfg)
rk <- merge(ao, cfg[, c(COMBO, "y"), with = FALSE], by = COMBO)
rk[, Truth := factor(V8_TRUTH_LABELS[truth], V8_TRUTH_LABELS[V8_KO_TRUTHS])]
print(cfg[1:8, .(label, best, worst, mean_rank = round(mean_rank, 2))])
write_sd(rk, SRC_DIR, STEM, "ranks")
write_sd(cfg, SRC_DIR, STEM, "configs")

# ---- BCFtools margin on each aorta truth set ----------------------------------------------
best <- ao[order(-mean_f1, overall_rank)][, .SD[1], by = .(truth, caller)]
m <- best[, {
  b <- .SD[caller == "bcftools"]; o <- .SD[caller != "bcftools"][order(-mean_f1)][1]
  .(bcftools_config = b$config, bcftools_f1 = b$mean_f1, other = o$caller, other_config = o$config,
    other_f1 = o$mean_f1, margin = b$mean_f1 - o$mean_f1)
}, by = truth]
SHORT <- setNames(c(sprintf("Power-gated, \u2265 %d%%", c(25, 50, 75)), sprintf("Depth-gated, \u2265 %d%%", c(25, 50, 75))),
                  V8_KO_TRUTHS)
m[, `:=`(Truth = factor(SHORT[truth], rev(SHORT)), Other = relabel(other, label_caller))]
print(m[, .(truth, bcftools_config, other_config, margin = round(margin, 4))])
stopifnot(m[, all(margin > 0)])
write_sd(m, SRC_DIR, STEM, "margin")

# ---- cross-species mean F1 ------------------------------------------------------------
cross_one <- function(d, ds) {
  r <- d[, .(mean_f1 = mean(mean_f1)), by = COMBO]
  r[, rank := frank(-mean_f1, ties.method = "min")][, dataset := ds]
}
cr <- dcast(rbind(cross_one(te, "teratoma"), cross_one(ao, "aorta")), ... ~ dataset, value.var = c("mean_f1", "rank"))
stopifnot(nrow(cr) == N_CFG)
add_labels(cr)
cr[, top6_both := rank_teratoma <= 6 & rank_aorta <= 6]
rho <- cor(cr$mean_f1_teratoma, cr$mean_f1_aorta, method = "spearman")
cat(sprintf("\nSpearman's rho, mean F1 over six KO truths, teratoma vs aorta: %.3f\n", rho))
print(cr[order(rank_teratoma)][1:8, .(label, rank_teratoma, rank_aorta)])
stopifnot(cr[rank_teratoma <= 6, all(rank_aorta <= 6 & caller == "bcftools" & strand_mode == "annotation")])
cr[, rho := rho]
write_sd(cr, SRC_DIR, STEM, "cross")
