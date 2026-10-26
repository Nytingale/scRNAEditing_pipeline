#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_rank_robustness_source_data.R   (2026-10-13; split from the _plot.R script 2026-10-25)
#
# Rank of every configuration under the nine v8 truth sets (definition: human/pm/gt_eval_v8/README.md;
# paths in plotting/gt_v8_paths.R):
#   call sets  nearby-mismatch pass only (*_step5_nomm; state mismatch)
#   KO truths  v8 (v7 rules, nearby-mismatch filter w1 = 200, PGP1 SNV, dbSNP common and
#              multiallelic sites removed, no region mask), expected / depth gate x t = 25 / 50 / 75%
#   REDIportal cov10_alt1 / alt3 / alt3_e05, PGP1 SNV, dbSNP common and multiallelic sites removed
# Rank = rank by mean F1 across arm / sample / cell-type tiers; ties broken by the sum of the per-tier F1 ranks.
#
# INPUT : human/pm/gt_eval_v8/compare_config_rank.tsv
#         human/pm/gt_eval_v8/scoring/mismatch/config_rank.tsv (parity)
#         human/pm/gt_eval_v8/rediportal/scoring_mismatch/gt_tier_scoring_*.tsv
# OUTPUT: source_data/pm_gt/TERATOMA_v8_rank_robustness.{ranks,configs,by_family}.tsv
#           ranks      rank of each configuration under each truth set
#           configs    configurations in figure order (worst rank, then mean rank) with best / worst / mean rank
#           by_family  the same within the six KO ground truths (ko) and the three REDIportal truths (redi)
#         (drawn by TERATOMA_v8_rank_robustness_plot.R: response figures 4 and 5)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_rank_robustness"

STATE <- V8_STATE
CMP <- V8_CMP_RANK
REF_RANK <- v8_config_rank(STATE)
REDI_V <- V8_REDI_VARIANTS
REDI_FILES <- unname(v8_redi_files(STATE))
stopifnot(all(file.exists(c(CMP, REF_RANK, REDI_FILES))))

COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
TRUTH_LAB <- V8_TRUTH_LABELS[V8_TRUTH_IDS]

cmp <- fread(CMP)
ref <- fread(REF_RANK)
chk <- merge(cmp[state == STATE, .(gate, t, config, a = overall_rank)],
             ref[, .(gate, t, config, b = overall_rank)], by = c("gate", "t", "config"))
stopifnot(nrow(chk) == nrow(ref), chk[, all(a == b)])
rk_ko <- cmp[state == STATE]
rk_ko[, truth := paste(gate, t, sep = "_")]
rk_ko <- rk_ko[, c(COMBO, "rank_sum", "mean_f1", "overall_rank", "truth"), with = FALSE]

overall_rank_redi <- function(s) {
  s <- copy(s)   # f1 = mean of the unit F1 (gt_score_rediportal_v8.R), as the KO truths
  s[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
  r <- s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
  r[, overall_rank := frank(list(-mean_f1, rank_sum), ties.method = "min")]
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
print(cfg_rb[1:15, .(label, best, worst, mean_rank = round(mean_rank, 1))])
write_sd(rk, SRC_DIR, STEM, "ranks")

N <- nrow(cfg_rb)
all_cfg <- copy(cfg_rb)[, y := N + 1L - .I]
all_cfg[, `:=`(Preprocessing = relabel(bam_type, label_bamType), Caller = relabel(caller, label_caller),
               Stranding = relabel(strand_mode, label_strandMode), `SNP removal` = relabel(snp_subtract, label_snpSubtract))]
write_sd(all_cfg, SRC_DIR, STEM, "configs")

# Within each truth family, configurations sorted by the worst rank within the family.
FAMILIES <- list(ko = V8_KO_TRUTHS, redi = paste0("redi_", REDI_V))
fam_rank <- rbindlist(lapply(names(FAMILIES), function(nm) {
  cfg <- rk[truth %in% FAMILIES[[nm]], .(worst = max(overall_rank), best = min(overall_rank),
                                         mean_rank = mean(overall_rank)), by = c(COMBO, "label")]
  setorder(cfg, worst, mean_rank)
  cfg[, `:=`(y = .N + 1L - .I, Preprocessing = relabel(bam_type, label_bamType), Caller = relabel(caller, label_caller),
             Stranding = relabel(strand_mode, label_strandMode), `SNP removal` = relabel(snp_subtract, label_snpSubtract))]
  cfg[, fam := nm]
}))
fam_w <- dcast(fam_rank, label ~ fam, value.var = "mean_rank")
cat("\nSpearman's rho, mean rank within KO ground truths vs within REDIportal:",
    round(cor(fam_w$ko, fam_w$redi, method = "spearman"), 3), "\n")
write_sd(fam_rank, SRC_DIR, STEM, "by_family")
