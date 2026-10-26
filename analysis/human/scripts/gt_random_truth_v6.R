#!/usr/bin/env Rscript
# =============================================================================
# gt_random_truth_v6.R   (2026-09-29)
#
# Random-truth control for the v6 scoring ("not halved = FP"). The labelled
# universe is the 8,553 v6 sites with a label (TRUE + FP_known_not_depleted +
# FP_novel). The 1,692 TRUE labels are replaced by random draws of the same
# size and every NTC call set is re-scored exactly as gt_score_v6.R does:
#   precision = TP / labelled calls, recall = TP / 1,692, F1 = 2PR / (P + R),
#   tier value = mean over units (arm 1, sample 4, cell type n).
# Two random designs, N_DRAW draws each:
#   matched  within strata of NTC editing x NTC coverage x KO information
#            (lambda_max = max over arms of KO cov x NTC editing), each stratum
#            gets as many random TRUE as the real set has there
#   permuted plain permutation of the labels over the universe
# If a caller ranking is reproduced by matched random truth, the ranking reflects
# where the calls sit (depth, editing level), not KO depletion.
#
# Inputs : gt_eval_v6/truth/gt_truth_sites.tsv.gz
#          gt_eval_v6/scoring/<tier>/gt_v6_calls_NTC.tsv.gz (from gt_score_v6.R)
# Outputs: gt_eval_v6/scoring/random_truth/
#          strata.tsv              stratum sizes and real TRUE counts
#          scores_long.tsv.gz      tier x config x truth draw: precision, recall, F1
#          config_summary.tsv      real vs random mean / 2.5% / 97.5% per config
#          rank_agreement.tsv      Spearman of config F1 / precision, real vs each draw
#          caller_best.tsv         best reported config per caller, real vs random
#          random_truth.pdf/.png
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(ggplot2); library(cowplot) })
B = Sys.getenv("SCRNA_DIR", ".")
V6 = file.path(B, "human/pm/gt_eval_v6")
OUT = file.path(V6, "scoring", "random_truth"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(B, "plotting/plotting_functions.R"))
N_DRAW = 200
TIERS = c("arm", "sample", "celltype")
LAB = c("TRUE", "FP_known_not_depleted", "FP_novel")
UNREPORTED_BAM = "sample"
# A BAM-processing configuration whose defining parameter has no effect on a caller's calls is
# recorded but not reported for that caller (HANDOFF 2026-09-29 UPDATE-19). scAllele has no per-base
# quality filter: recalq20 == recal, and BQSR changes only FILTER, so its recal BAM differs from
# sampleRmDup only by SplitNCigarReads.
UNREPORTED_CALLER_BAM = list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))
unreported = function(caller, bam_type) bam_type %in% UNREPORTED_BAM |
  mapply(function(c, b) b %in% UNREPORTED_CALLER_BAM[[c]], caller, bam_type) &
  caller %in% names(UNREPORTED_CALLER_BAM)
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")
set.seed(20260929)

# ---- universe and strata ---------------------------------------------------------
tr = fread(file.path(V6, "truth", "gt_truth_sites.tsv.gz"),
           select = c("site", "class", "cov_NTC_ctrl", "edit_NTC_ctrl", "cov_ADAR1_KO", "cov_ADARB1_KO"))
U = tr[class %in% LAB]
U[, lambda_max := pmax(cov_ADAR1_KO, cov_ADARB1_KO) * edit_NTC_ctrl]
U[, `:=`(s_edit = cut(edit_NTC_ctrl, c(0, 0.1, 0.2, 0.35, 0.5, 1.01), right = FALSE),
         s_cov = cut(cov_NTC_ctrl, unique(quantile(cov_NTC_ctrl, 0:4 / 4)), include.lowest = TRUE),
         s_lam = cut(lambda_max, c(0, 3, 5, 10, 20, Inf), right = FALSE))]
U[, stratum := paste(s_edit, s_cov, s_lam, sep = "|")]
U[, idx := .I]
is_true = U$class == "TRUE"
N_TRUE = sum(is_true)
st = U[, .(n = .N, n_true = sum(class == "TRUE")), by = .(stratum, s_edit, s_cov, s_lam)][order(stratum)]
fwrite(st, file.path(OUT, "strata.tsv"), sep = "\t")
msg("universe ", nrow(U), " sites, ", N_TRUE, " TRUE, ", nrow(st), " strata; ",
    st[n_true == n & n > 0, sum(n)], " TRUE sites in strata that are all TRUE (fixed in every matched draw)")

by_str = split(U$idx, U$stratum); k_str = st[match(names(by_str), stratum), n_true]
draw_matched = function() { v = logical(nrow(U)); for (j in seq_along(by_str)) if (k_str[j] > 0)
  v[by_str[[j]][sample.int(length(by_str[[j]]), k_str[j])]] = TRUE; v }
Tm = cbind(real = is_true,
           sapply(seq_len(N_DRAW), function(i) draw_matched()),
           sapply(seq_len(N_DRAW), function(i) sample(is_true)))
colnames(Tm) = c("real", sprintf("matched_%03d", seq_len(N_DRAW)), sprintf("permuted_%03d", seq_len(N_DRAW)))
stopifnot(all(colSums(Tm) == N_TRUE))
ov = colSums(Tm[, -1] & is_true) / N_TRUE
msg("overlap with real TRUE: matched mean ", round(mean(ov[grep("^matched", names(ov))]), 3),
    ", permuted mean ", round(mean(ov[grep("^permuted", names(ov))]), 3))
Tm = Matrix(Tm * 1, sparse = FALSE)

# ---- re-score every tier ---------------------------------------------------------
long = rbindlist(lapply(TIERS, function(tier) {
  cl = fread(file.path(V6, "scoring", tier, "gt_v6_calls_NTC.tsv.gz"))
  us = fread(file.path(V6, "scoring", tier, "gt_v6_unit_scores.tsv"), select = c("unit", "config", "caller", "bam_type"))
  us = us[!unreported(caller, bam_type)]
  cl = cl[config %in% us$config & class %in% LAB]
  cl[, idx := U$idx[match(site, U$site)]]
  stopifnot(!anyNA(cl$idx))
  us[, row := .I]
  cl[, row := us$row[match(paste(unit, config), paste(us$unit, us$config))]]
  cl = cl[!is.na(row)]
  M = sparseMatrix(i = cl$row, j = cl$idx, x = 1, dims = c(nrow(us), nrow(U)))
  tp = as.matrix(M %*% Tm); nl = rowSums(M)
  P = tp / nl; P[nl == 0, ] = NA; R = tp / N_TRUE
  F1 = ifelse(!is.na(P) & P + R > 0, 2 * P * R / (P + R), NA)
  d = data.table(unit = rep(us$unit, ncol(Tm)), config = rep(us$config, ncol(Tm)), caller = rep(us$caller, ncol(Tm)),
                 truth = rep(colnames(Tm), each = nrow(us)),
                 precision = as.vector(P), recall = as.vector(R), f1 = as.vector(F1))
  msg(tier, ": ", nrow(us), " unit x config rows, ", uniqueN(us$config), " reported configs")
  d[, .(precision = mean(precision, na.rm = TRUE), recall = mean(recall), f1 = mean(f1, na.rm = TRUE)),
    by = .(config, caller, truth)][, tier := tier]
}))
long[, design := sub("_[0-9]+$", "", truth)]
fwrite(long, file.path(OUT, "scores_long.tsv.gz"), sep = "\t")

# real scores must equal gt_score_v6.R's reported scores
chk = rbindlist(lapply(TIERS, function(tier)
  fread(file.path(V6, "scoring", tier, "gt_v6_scores_reported.tsv"), select = c("config", "f1"))[, tier := tier]))
chk = merge(chk, long[truth == "real", .(tier, config, f1_here = f1)], by = c("tier", "config"))
msg("parity with gt_score_v6.R: max |dF1| = ", signif(max(abs(chk$f1 - chk$f1_here), na.rm = TRUE), 3), " over ", nrow(chk), " configs")

# ---- summaries -------------------------------------------------------------------
real = long[truth == "real", .(tier, config, caller, P_real = precision, R_real = recall, F_real = f1)]
rnd = long[design != "real", .(P_mean = mean(precision, na.rm = TRUE), P_lo = quantile(precision, 0.025, na.rm = TRUE),
                                P_hi = quantile(precision, 0.975, na.rm = TRUE),
                                R_mean = mean(recall), F_mean = mean(f1, na.rm = TRUE),
                                F_lo = quantile(f1, 0.025, na.rm = TRUE), F_hi = quantile(f1, 0.975, na.rm = TRUE),
                                F_sd = sd(f1, na.rm = TRUE)),
           by = .(tier, design, config)]
cs = merge(real, rnd, by = c("tier", "config"))
cs[, `:=`(F_excess = F_real - F_mean, F_z = (F_real - F_mean) / F_sd)]
setorder(cs, tier, design, -F_real)
fwrite(cs, file.path(OUT, "config_summary.tsv"), sep = "\t")

rk = long[design != "real", {
  cfg = config
  r = real[tier == .BY$tier]; r = r[match(cfg, r$config)]
  .(rho_f1 = cor(f1, r$F_real, method = "spearman", use = "complete.obs"),
    rho_precision = cor(precision, r$P_real, method = "spearman", use = "complete.obs"),
    rho_recall = cor(recall, r$R_real, method = "spearman", use = "complete.obs"))
}, by = .(tier, design, truth)]
fwrite(rk, file.path(OUT, "rank_agreement.tsv"), sep = "\t")

best = long[, .SD[which.max(f1)], by = .(tier, truth, caller)]
best[, rank := frank(-f1, ties.method = "min"), by = .(tier, truth)]
real_order = best[truth == "real", .(tier, caller, rank_real = rank)]
cb = merge(best[design != "real"], real_order, by = c("tier", "caller"))
same_order = cb[, .(same = all(rank == rank_real)), by = .(tier, design, truth)][, .(frac_draws_same_caller_order = mean(same)), by = .(tier, design)]
cbs = merge(best[truth == "real", .(tier, caller, config_real = config, F_real = f1, rank_real = rank)],
            cb[, .(F_rand_mean = mean(f1), F_rand_lo = quantile(f1, 0.025), F_rand_hi = quantile(f1, 0.975),
                   rank_rand_mean = mean(rank), p_same_rank = mean(rank == rank_real)), by = .(tier, design, caller)],
            by = c("tier", "caller"))
setorder(cbs, tier, design, rank_real)
fwrite(cbs, file.path(OUT, "caller_best.tsv"), sep = "\t")
print(same_order)
print(rk[, .(rho_f1_median = median(rho_f1), rho_f1_lo = quantile(rho_f1, 0.025), rho_prec_median = median(rho_precision)),
         by = .(tier, design)])
print(cbs[, .(tier, design, caller, F_real = round(F_real, 3), F_rand = round(F_rand_mean, 3),
              rank_real, rank_rand = round(rank_rand_mean, 2), p_same_rank = round(p_same_rank, 2))])

# ---- figure ----------------------------------------------------------------------
cs[, `:=`(Caller = std_caller(caller), Tier = factor(tier, TIERS, c("Arm", "Sample", "Cell type")),
          Design = factor(design, c("matched", "permuted"), c("Matched random truth", "Permuted random truth")))]
cols = cmap_method[unique(cs$Caller)]
pa = ggplot(cs, aes(F_mean, F_real, colour = Caller)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = F_lo, xmax = F_hi), height = 0, alpha = 0.4) +
  geom_point(size = 1.8) + scale_colour_manual(values = cols) +
  facet_grid(Design ~ Tier, scales = "free") +
  labs(x = "F1 with random TRUE labels (mean, 95% of draws)", y = "F1 with real v6 TRUE labels", colour = NULL) +
  theme_bw(11) + theme(legend.position = "bottom")
pp = ggplot(cs, aes(P_mean, P_real, colour = Caller)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = P_lo, xmax = P_hi), height = 0, alpha = 0.4) +
  geom_point(size = 1.8) + scale_colour_manual(values = cols) +
  facet_grid(Design ~ Tier, scales = "free") +
  labs(x = "Precision with random TRUE labels", y = "Precision with real v6 TRUE labels", colour = NULL) +
  theme_bw(11) + theme(legend.position = "none")
rk[, `:=`(Tier = factor(tier, TIERS, c("Arm", "Sample", "Cell type")), Design = factor(design, c("matched", "permuted"), c("Matched", "Permuted")))]
pr = ggplot(rk, aes(Design, rho_f1)) + geom_boxplot(outlier.size = 0.5) + facet_wrap(~Tier, nrow = 1) +
  labs(x = NULL, y = "Spearman rho of config F1\n(real vs one random draw)") + theme_bw(11)
p = plot_grid(pa, plot_grid(pp, pr, nrow = 1, rel_widths = c(1.4, 1), labels = c("b", "c")), ncol = 1,
              labels = c("a", ""), rel_heights = c(1.15, 1))
ggsave(file.path(OUT, "random_truth.pdf"), p, width = 15, height = 13)
ggsave(file.path(OUT, "random_truth.png"), p, width = 15, height = 13, dpi = 90)
msg("done -> ", OUT)
