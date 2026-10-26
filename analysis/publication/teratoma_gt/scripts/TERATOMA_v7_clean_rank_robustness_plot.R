#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v7_clean_rank_robustness_plot.R   (2026-10-11)
#
# TERATOMA_v7_rank_robustness_all_tiles (superseded 2026-10-12, ../_superseded/v7_uncleaned_rank_2026-10-12/)
# with cleaned truths and
# the call sets without the "retain clustered sites" filter:
#   call sets  nearby-mismatch pass only (gt_eval_v6/prewindow/results_<tier> *_step5_nomm)
#   KO truths  v7_mask_nomm_pgp1 (gt_truth_v7_variants.R): artifact-prone mask + truth-level
#              nearby-mismatch filter, PGP1 SNV positions removed; multiallelic and dbSNP common
#              sites removed as in v7; expected / depth gate x t = 25 / 50 / 75%
#   REDIportal cov10_alt1 / alt3 / alt3_e05 (gt_score_rediportal_v6.R): PGP1 SNV positions,
#              dbSNP common and multiallelic sites removed from the positives and the candidate
#              universe, then the three gates (v7 default; run_gt_rediportal_v7_states.sh mismatch)
#
# INPUT : human/pm/gt_eval_v7/variants/compare_config_rank.tsv (gt_v7_variants_compare.R)
#         human/pm/gt_eval_v7/scoring/config_rank.tsv (parity of the rank computation)
#         human/pm/gt_eval_v7/rediportal/scoring_mismatch/gt_tier_scoring_*.tsv
# OUTPUT: figures/<DATASET>/TERATOMA_v7_clean_rank_robustness_all_tiles{,_legend}.svg
#         source_data/<DATASET>/TERATOMA_v7_clean_rank_robustness.tsv
# Usage : Rscript TERATOMA_v7_clean_rank_robustness_plot.R [DATASET]
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "TERATOMA_common.R"))
argv <- commandArgs(trailingOnly = TRUE)
teratoma_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

VARIANT <- "v7_mask_nomm_pgp1"
STATE <- "mismatch"
V7 <- file.path(REPO_DIR, "human", "pm", "gt_eval_v7")
CMP <- file.path(V7, "variants", "compare_config_rank.tsv")
REF_RANK <- file.path(V7, "scoring", "config_rank.tsv")
REDI_DIR <- file.path(V7, "rediportal", "scoring_mismatch")
REDI_V <- c("cov10_alt1", "alt3", "alt3_e05")
REDI_FILES <- file.path(REDI_DIR, paste0("gt_tier_scoring_", REDI_V, ".tsv"))
SRC <- c(CMP, REF_RANK, REDI_FILES)
stopifnot(all(file.exists(SRC)))

COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
V7_TRUTHS <- data.table(
  truth = c("expected_t25", "expected_t50", "expected_t75", "depth_t25", "depth_t50", "depth_t75",
            paste0("redi_", REDI_V)),
  label = NINE_TRUTHS$label, colour = NINE_TRUTHS$colour)
TRUTH_LAB <- setNames(V7_TRUTHS$label, V7_TRUTHS$truth)
cmap_truth <- setNames(V7_TRUTHS$colour, V7_TRUTHS$label)
N_TRUTH <- nrow(V7_TRUTHS)

# ---- ranks ---------------------------------------------------------------------
cmp <- fread(CMP)
ref <- fread(REF_RANK)
chk <- merge(cmp[variant == "v7" & state == "both", .(gate, t, config, a = overall_rank)],
             ref[, .(gate, t, config, b = overall_rank)], by = c("gate", "t", "config"))
stopifnot(nrow(chk) == nrow(ref), chk[, all(a == b)])
rk_ko <- cmp[variant == VARIANT & state == STATE]
rk_ko[, truth := paste(gate, t, sep = "_")]
rk_ko <- rk_ko[, c(COMBO, "rank_sum", "mean_f1", "overall_rank", "truth"), with = FALSE]

overall_rank_redi <- function(s) {
  s <- copy(s)[, f1 := fifelse(precision + recall_vs_union > 0,
                               2 * precision * recall_vs_union / (precision + recall_vs_union), NA_real_)]
  s[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
  r <- s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
  r[, overall_rank := frank(list(rank_sum, -mean_f1), ties.method = "min")]
}
rk_redi <- rbindlist(lapply(seq_along(REDI_V), function(i)
  overall_rank_redi(fread(REDI_FILES[i]))[, truth := paste0("redi_", REDI_V[i])]))
rk <- rbind(rk_ko, rk_redi)
stopifnot(rk[, .N, by = truth][, all(N == N[1])], uniqueN(rk[, ..COMBO]) == rk[, .N, by = truth]$N[1])
rk[, Truth := factor(TRUTH_LAB[truth], TRUTH_LAB)]
cfg_rb <- rk[, .(worst = max(overall_rank), best = min(overall_rank), mean_rank = mean(overall_rank),
                 mean_rank_sum = mean(rank_sum)), by = COMBO]
setorder(cfg_rb, worst, mean_rank)
cfg_rb[, label := sprintf("%s | %s | %s | %s", relabel(caller, label_caller), relabel(bam_type, label_bamType),
                          relabel(strand_mode, label_strandMode), relabel(snp_subtract, label_snpSubtract))]
cfg_rb[, label := factor(label, rev(label))]
rk <- merge(rk, cfg_rb, by = COMBO)
fwrite(rk[order(worst, mean_rank, Truth)], file.path(SRC_DIR, "TERATOMA_v7_clean_rank_robustness.tsv"), sep = "\t")
print(cfg_rb[1:15, .(label, best, worst, mean_rank = round(mean_rank, 1))])

# ---- tiles figure (as TERATOMA_v7_rank_robustness_all_tiles) --------------------
N <- nrow(cfg_rb)
all_cfg <- copy(cfg_rb)[, y := N + 1L - .I]
all_cfg[, `:=`(Preprocessing = relabel(bam_type, label_bamType), Caller = relabel(caller, label_caller),
               Stranding = relabel(strand_mode, label_strandMode), `SNP removal` = relabel(snp_subtract, label_snpSubtract))]
CMAPS <- list(Preprocessing = cmap_proc, Caller = cmap_method, Stranding = cmap_strand, `SNP removal` = cmap_snp)
Y_LIM_ALL <- c(0.5, N + 0.5)
BASE_ALL <- 20
tile_col <- function(col) {
  d <- all_cfg[, .(y, value = get(col), column = col)]
  ggplot(d, aes(column, y, fill = value)) +
    geom_tile(width = 1, height = 0.78) +
    scale_fill_manual(values = CMAPS[[col]], name = col,
                      breaks = if (col == "Caller") {
                        v <- unique(d$value); v[order(tolower(v))]
                      } else intersect(names(CMAPS[[col]]), unique(d$value))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 1, title.position = "top")) +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = BASE_ALL) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
          axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(5, 0, 5, 0))
}
d <- copy(rk)
d[, y := (N + 1L - match(as.character(label), as.character(cfg_rb$label))) +
        (as.integer(Truth) - (N_TRUTH + 1) / 2) * 0.04]
cf <- copy(all_cfg)
p_dots_all <- ggplot(d, aes(overall_rank, y)) +
  geom_segment(data = cf, aes(x = best, xend = worst, y = y, yend = y), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = Truth), size = 4.6) +
  geom_text(data = cf, aes(x = worst, y = y, label = worst), nudge_x = 0.9, hjust = 0, size = 7, fontface = "bold") +
  scale_x_continuous(breaks = c(1, seq(5, 55, 5)), expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_continuous(limits = Y_LIM_ALL, expand = c(0, 0)) +
  scale_colour_manual_cvd(values = cmap_truth, name = NULL) +
  guides(colour = guide_legend(ncol = 1, title = "Truth set", title.position = "top", override.aes = list(size = 5))) +
  labs(x = "Overall cross-tier rank", y = NULL) +
  theme_manuscript(base_size = BASE_ALL) + grid_x +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey82"), plot.margin = margin(5, 5, 5, 4))
TILE_IN <- 0.36
legend_theme <- function(text_size, key_cm, just = "top")
  theme(legend.position = "bottom", legend.justification = just,
        legend.title = element_text(face = "bold", size = text_size + 1),
        legend.text = element_text(size = text_size), legend.key.size = unit(key_cm, "cm"))
panels <- c(lapply(names(CMAPS), tile_col), list(p_dots_all))
p_tiles <- patchwork::wrap_plots(panels, nrow = 1,
  widths = unit(c(rep(TILE_IN, length(CMAPS)), 1), c(rep("in", length(CMAPS)), "null"))) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = paste0("Rank of ", N, " method configurations under ", N_TRUTH, " truth definitions"),
    subtitle = paste0("Call sets: nearby-mismatch pass only (no clustered-sites filter). KO truths: mask + nearby-mismatch ",
                      "filter, PGP1 SNVs removed.\nREDIportal: PGP1 SNVs and multiallelic sites removed. All truths: ",
                      "dbSNP common and multiallelic sites removed"),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
                  plot.subtitle = element_text(hjust = 0.5, size = 14))) &
  legend_theme(10, 0.35)
save_figure(p_tiles, "TERATOMA_v7_clean_rank_robustness_all_tiles", width = 16, height = 3 + TILE_IN * N + 2.2,
            sources = SRC, data = all_cfg, grid = "x", script = SCRIPT)

legs <- lapply(seq_along(panels), function(i) {
  pp <- panels[[i]] + legend_theme(18, 0.8, just = "left")
  if (i == length(panels))
    pp <- pp + guides(colour = guide_legend(ncol = 1, title = "Truth set", title.position = "top",
                                            override.aes = list(size = 7)))
  cowplot::get_plot_component(pp, "guide-box-bottom")
})
n_keys <- c(sapply(names(CMAPS), function(col) uniqueN(all_cfg[[col]])), N_TRUTH)
leg_h <- n_keys + 2.2
p_leg <- cowplot::plot_grid(plotlist = legs, ncol = 1, rel_heights = leg_h, align = "v", axis = "l")
save_figure(p_leg, "TERATOMA_v7_clean_rank_robustness_all_tiles_legend", width = 8, height = 0.4 * sum(leg_h),
            sources = SRC, data = all_cfg, grid = "none", script = SCRIPT)
message("Done")
