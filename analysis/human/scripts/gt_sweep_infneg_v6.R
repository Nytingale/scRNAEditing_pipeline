#!/usr/bin/env Rscript
# =============================================================================
# gt_sweep_infneg_v6.R   (2026-09-29)
#
# Depletion-threshold sweep scored against INFORMATIVE negatives, with a
# random-truth control (HANDOFF_2026-09-29 UPDATE-16/17).
#
# Eligible universe E: the v6 gates (evaluable, NTC edited, reproducible) = 8,553 sites.
# Per KO arm (ADAR1_KO, ADARB1_KO), with the gt_eval_classify_v6.R definitions:
#   covered_a  NTC cov >= 10 and KO cov x NTC editing >= 3
#   dep_a(t)   covered_a and KO editing <= (1 - t) x NTC editing
# Classes at threshold t in {0.25, 0.50, 0.75}:
#   positive               E and dep_a(t) in either arm (t = 0.50 is v6 TRUE)
#   informative negative   E, covered in BOTH arms, and the one-sided 95%
#                          Clopper-Pearson lower bound on KO editing > 0.75 x NTC
#                          editing in BOTH arms (confidently < 25% depleted);
#                          fixed across t; NTC editing treated as known
#   unresolved             everything else (not scored)
# Scoring per unit x configuration (A>G NTC calls, reported configurations):
#   TP = calls on positives, FP = calls on informative negatives,
#   precision = TP / (TP + FP), recall = TP / positives, F1; tier value = mean over units.
# Label-free call descriptives: fraction of calls scorable (positive or negative),
#   fraction not in the candidate table, fraction on candidates with >= 20 NTC
#   reads that fail the NTC-edited gate (weak NTC support).
# Random-truth designs, N_DRAW draws each, per t (strata: NTC editing x NTC
# coverage quartile x lambda_max, as gt_random_truth_v6.R):
#   resolved_matched   labels shuffled within positives + negatives, within strata
#   resolved_permuted  labels shuffled within positives + negatives
#   universe_matched   random positives and negatives drawn from E, disjoint,
#                      each matched to the real class's strata
#
# Inputs : gt_eval_v6/truth/gt_truth_sites.tsv.gz
#          gt_eval_v6/scoring/<tier>/gt_v6_calls_NTC.tsv.gz, gt_v6_unit_scores.tsv
# Outputs: gt_eval_v6/scoring/sweep_infneg/
#          truth_sweep_classes.tsv.gz (per eligible site: class at each t),
#          unit_scores_real.tsv.gz (per unit x config x t, real labels),
#          arm_bootstrap_se.tsv (arm tier: site-bootstrap SE over resolved sites, N_BOOT reps),
#          class_counts.tsv, scores_long.tsv.gz, call_descriptives.tsv,
#          config_summary.tsv, rank_agreement.tsv, caller_best.tsv,
#          sweep_rank_stability.tsv, sweep_infneg.{pdf,png}
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(ggplot2); library(cowplot) })
B = Sys.getenv("SCRNA_DIR", ".")
V6 = file.path(B, "human/pm/gt_eval_v6")
# GT_SCORING_DIR: scoring root holding <tier>/gt_v6_calls_NTC.tsv.gz (default gt_eval_v6/scoring;
# the pre-window comparison uses gt_eval_v6/prewindow/scoring)
SCORING = Sys.getenv("GT_SCORING_DIR", file.path(V6, "scoring"))
# GT_COVER_GATE: how a KO arm counts as covered (nine-truth suite, 2026-10-01)
#   expected (default, frozen v6)  NTC cov >= 10 and KO cov x NTC editing >= 3 -> sweep_infneg/
#   depth                          NTC cov >= 10 and KO cov >= 10              -> sweep_infneg_depth/
#   E is rebuilt from the gate: ctrl_edited and consistent and covered in either arm.
COVER_GATE = Sys.getenv("GT_COVER_GATE", "expected")
stopifnot(COVER_GATE %in% c("expected", "depth"))
OUT = file.path(SCORING, if (COVER_GATE == "expected") "sweep_infneg" else "sweep_infneg_depth")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(B, "plotting/plotting_functions.R"))
N_DRAW = 200; N_BOOT = 500
TS = c(0.25, 0.50, 0.75)
NEG_MAX_DEP = 0.25; CI_ALPHA = 0.05; MIN_COV = 10; KO_MIN_EXP = 3; MIN_COV_KO_DEPTH = 10; WEAK_MIN_COV = 20
TIERS = c("arm", "sample", "celltype")
ARMS = c("ADAR1_KO", "ADARB1_KO")
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

# ---- classes ---------------------------------------------------------------------
tr = fread(file.path(V6, "truth", "gt_truth_sites.tsv.gz"))
p = tr$edit_NTC_ctrl
cp_lower = function(k, n, a = CI_ALPHA) { o = rep(0, length(k)); ok = n > 0 & k > 0; o[ok] = qbeta(a, k[ok], n[ok] - k[ok] + 1); o }
covered = sapply(ARMS, function(a) tr$cov_NTC_ctrl >= MIN_COV &
  if (COVER_GATE == "expected") tr[[paste0("cov_", a)]] * p >= KO_MIN_EXP else tr[[paste0("cov_", a)]] >= MIN_COV_KO_DEPTH)
E = (covered[, 1] | covered[, 2]) & tr$ctrl_edited & tr$consistent
if (COVER_GATE == "expected") stopifnot(identical(E, tr$evaluable & tr$ctrl_edited & tr$consistent))
neg = E & covered[, 1] & covered[, 2] &
  cp_lower(tr$alt_ADAR1_KO, tr$cov_ADAR1_KO) > (1 - NEG_MAX_DEP) * p &
  cp_lower(tr$alt_ADARB1_KO, tr$cov_ADARB1_KO) > (1 - NEG_MAX_DEP) * p
pos = sapply(TS, function(t) E & ((covered[, 1] & tr$edit_ADAR1_KO <= (1 - t) * p) |
                                  (covered[, 2] & tr$edit_ADARB1_KO <= (1 - t) * p)))
colnames(pos) = sprintf("t%02d", 100 * TS)
stopifnot(!any(pos & neg))
if (COVER_GATE == "expected") stopifnot(sum(pos[, "t50"]) == sum(tr$class == "TRUE"))
msg("cover gate: ", COVER_GATE, "; eligible ", sum(E), "; output ", OUT)
weak = !E & tr$cov_NTC_ctrl >= WEAK_MIN_COV & !tr$ctrl_edited
cc = data.table(t = TS, positive = colSums(pos), informative_negative = sum(neg), eligible = sum(E),
                unresolved = sum(E) - colSums(pos) - sum(neg))
fwrite(cc, file.path(OUT, "class_counts.tsv"), sep = "\t"); print(cc)
cls_tab = data.table(site = tr$site, eligible = E, covered_ADAR1 = covered[, 1], covered_ADARB1 = covered[, 2],
                     informative_negative = neg)
for (j in colnames(pos)) set(cls_tab, j = paste0("class_", j),
  value = fifelse(!E, "not_eligible", fifelse(pos[, j], "positive", fifelse(neg, "informative_negative", "unresolved"))))
fwrite(cls_tab[eligible == TRUE], file.path(OUT, "truth_sweep_classes.tsv.gz"), sep = "\t")

# ---- strata over E -----------------------------------------------------------------
U = tr[E, .(site, cov_NTC_ctrl, edit_NTC_ctrl, lambda_max = pmax(cov_ADAR1_KO, cov_ADARB1_KO) * edit_NTC_ctrl)]
U[, stratum := paste(cut(edit_NTC_ctrl, c(0, 0.1, 0.2, 0.35, 0.5, 1.01), right = FALSE),
                     cut(cov_NTC_ctrl, unique(quantile(cov_NTC_ctrl, 0:4 / 4)), include.lowest = TRUE),
                     cut(lambda_max, c(0, 3, 5, 10, 20, Inf), right = FALSE), sep = "|")]
posE = pos[E, , drop = FALSE]; negE = neg[E]
by_str = split(seq_len(nrow(U)), U$stratum)

shuffle_within = function(lab, set, groups) {    # permute lab among `set`, within groups
  out = lab
  for (g in split(set, groups[set])) if (length(g) > 1) out[g] = lab[g][sample.int(length(g))]
  out
}
draw_universe = function(isP, isN) {              # disjoint random P, N matched to real strata
  P = N = logical(length(isP))
  for (g in by_str) {
    kP = sum(isP[g]); kN = sum(isN[g])
    if (kP + kN == 0) next
    s = g[sample.int(length(g), kP + kN)]
    if (kP) P[s[seq_len(kP)]] = TRUE
    if (kN) N[s[kP + seq_len(kN)]] = TRUE
  }
  list(P, N)
}
labels = lapply(seq_along(TS), function(j) {
  isP = posE[, j]; res = which(isP | negE)
  cls = ifelse(isP, 1L, ifelse(negE, -1L, 0L))
  mk = function(design, f) lapply(seq_len(N_DRAW), function(i) f())
  d_rm = mk("resolved_matched", function() shuffle_within(cls, res, U$stratum))
  d_rp = mk("resolved_permuted", function() shuffle_within(cls, res, rep("all", nrow(U))))
  d_um = mk("universe_matched", function() { x = draw_universe(isP, negE); ifelse(x[[1]], 1L, ifelse(x[[2]], -1L, 0L)) })
  L = cbind(real = cls, do.call(cbind, d_rm), do.call(cbind, d_rp), do.call(cbind, d_um))
  colnames(L) = c("real", sprintf("resolved_matched_%03d", 1:N_DRAW), sprintf("resolved_permuted_%03d", 1:N_DRAW),
                  sprintf("universe_matched_%03d", 1:N_DRAW))
  stopifnot(all(colSums(L == 1L) == sum(isP)), all(colSums(L == -1L) == sum(negE)))
  L
})
names(labels) = colnames(pos)
msg("labels built: ", ncol(labels[[1]]), " truth columns per threshold")

# ---- score -------------------------------------------------------------------------
long = list(); desc = list(); unit_real = list(); arm_boot = list()
for (tier in TIERS) {
  cl = fread(file.path(SCORING, tier, "gt_v6_calls_NTC.tsv.gz"), select = c("unit", "config", "site"))
  us = unique(fread(file.path(SCORING, tier, "gt_v6_unit_scores.tsv"), select = c("unit", "config", "caller", "bam_type")))
  us = us[!unreported(caller, bam_type)][, row := .I]
  cl[, row := us$row[match(paste(unit, config), paste(us$unit, us$config))]]
  cl = cl[!is.na(row)]
  ti = match(cl$site, tr$site)
  n_calls = tabulate(cl$row, nrow(us))
  desc[[tier]] = data.table(tier = tier, us[, .(unit, config, caller)], n_calls = n_calls,
    frac_not_in_table = tabulate(cl$row[is.na(ti)], nrow(us)) / pmax(n_calls, 1),
    frac_weak_ntc = tabulate(cl$row[which(weak[ti])], nrow(us)) / pmax(n_calls, 1))
  ei = match(cl$site, U$site); k = !is.na(ei)
  M = sparseMatrix(i = cl$row[k], j = ei[k], x = 1, dims = c(nrow(us), nrow(U)))
  for (tn in names(labels)) {
    L = labels[[tn]]
    TP = as.matrix(M %*% Matrix(1 * (L == 1L))); FP = as.matrix(M %*% Matrix(1 * (L == -1L)))
    nP = colSums(L == 1L)[1]
    P = TP / (TP + FP); P[TP + FP == 0] = NA
    R = TP / nP
    F1 = ifelse(!is.na(P) & P + R > 0, 2 * P * R / (P + R), NA)
    d = data.table(unit = rep(us$unit, ncol(L)), config = rep(us$config, ncol(L)), caller = rep(us$caller, ncol(L)),
                   truth = rep(colnames(L), each = nrow(us)), tp = as.vector(TP), fp = as.vector(FP),
                   precision = as.vector(P), recall = as.vector(R), f1 = as.vector(F1))
    long[[paste(tier, tn)]] = d[, .(tp = mean(tp), fp = mean(fp), precision = mean(precision, na.rm = TRUE),
                                    recall = mean(recall), f1 = mean(f1, na.rm = TRUE)),
                                by = .(config, caller, truth)][, `:=`(tier = tier, t = tn)]
    unit_real[[paste(tier, tn)]] = d[truth == "real", !"truth"][, `:=`(tier = tier, t = tn, n_calls = n_calls)]
    if (tier == "arm") {   # site bootstrap over resolved sites, as gt_score_v6.R
      lab = L[, "real"]; res = which(lab != 0L)
      Mr = M[, res, drop = FALSE]; isP = lab[res] == 1L
      bs = replicate(N_BOOT, {
        i = sample.int(length(res), replace = TRUE)
        h = Mr[, i, drop = FALSE]
        tp = as.vector(h %*% as.numeric(isP[i])); fp = as.vector(h %*% as.numeric(!isP[i]))
        p = ifelse(tp + fp > 0, tp / (tp + fp), NA); r = tp / sum(isP[i])
        cbind(p, r, ifelse(!is.na(p) & p + r > 0, 2 * p * r / (p + r), NA))
      })
      arm_boot[[tn]] = data.table(tier = tier, t = tn, unit = us$unit, config = us$config,
        precision_err = apply(bs[, 1, ], 1, sd, na.rm = TRUE), recall_err = apply(bs[, 2, ], 1, sd, na.rm = TRUE),
        f1_err = apply(bs[, 3, ], 1, sd, na.rm = TRUE))
    }
  }
  msg(tier, ": ", nrow(us), " unit x config rows, ", uniqueN(us$config), " configs")
}
unit_real = rbindlist(unit_real)
for (v in c("precision", "f1")) set(unit_real, which(is.nan(unit_real[[v]])), v, NA_real_)
fwrite(unit_real, file.path(OUT, "unit_scores_real.tsv.gz"), sep = "\t")
fwrite(rbindlist(arm_boot), file.path(OUT, "arm_bootstrap_se.tsv"), sep = "\t")
long = rbindlist(long); long[, design := sub("_[0-9]+$", "", truth)]
for (v in c("precision", "f1")) set(long, which(is.nan(long[[v]])), v, NA_real_)
fwrite(long, file.path(OUT, "scores_long.tsv.gz"), sep = "\t")
desc = rbindlist(desc)[, .(n_calls = mean(n_calls), frac_not_in_table = mean(frac_not_in_table),
                           frac_weak_ntc = mean(frac_weak_ntc)), by = .(tier, config, caller)]
sc = long[truth == "real", .(tier, t, config, tp, fp)]
desc = merge(desc, dcast(sc, tier + config ~ t, value.var = c("tp", "fp")), by = c("tier", "config"))
fwrite(desc, file.path(OUT, "call_descriptives.tsv"), sep = "\t")

# ---- summaries -----------------------------------------------------------------------
real = long[truth == "real", .(tier, t, config, caller, P_real = precision, R_real = recall, F_real = f1, fp_real = fp)]
rnd = long[design != "real", .(P_mean = mean(precision, na.rm = TRUE), P_lo = quantile(precision, 0.025, na.rm = TRUE),
                                P_hi = quantile(precision, 0.975, na.rm = TRUE), R_mean = mean(recall),
                                F_mean = mean(f1, na.rm = TRUE), F_lo = quantile(f1, 0.025, na.rm = TRUE),
                                F_hi = quantile(f1, 0.975, na.rm = TRUE)), by = .(tier, t, design, config)]
cs = merge(real, rnd, by = c("tier", "t", "config"))
setorder(cs, tier, t, design, -F_real)
fwrite(cs, file.path(OUT, "config_summary.tsv"), sep = "\t")

rk = long[design != "real", {
  cfg = config; r = real[tier == .BY$tier & t == .BY$t]; r = r[match(cfg, r$config)]
  .(rho_f1 = cor(f1, r$F_real, method = "spearman", use = "complete.obs"),
    rho_precision = cor(precision, r$P_real, method = "spearman", use = "complete.obs"),
    rho_recall = cor(recall, r$R_real, method = "spearman", use = "complete.obs"))
}, by = .(tier, t, design, truth)]
fwrite(rk, file.path(OUT, "rank_agreement.tsv"), sep = "\t")

best = unique(long[!is.na(f1)][order(-f1)], by = c("tier", "t", "truth", "caller"))
best[, rank := frank(-f1, ties.method = "min"), by = .(tier, t, truth)]
ro = best[truth == "real", .(tier, t, caller, config_real = config, F_real = f1, P_real = precision, R_real = recall, rank_real = rank)]
cb = merge(best[design != "real"], ro[, .(tier, t, caller, rank_real)], by = c("tier", "t", "caller"))
cbs = merge(ro, cb[, .(F_rand_mean = mean(f1), P_rand_mean = mean(precision, na.rm = TRUE),
                       rank_rand_mean = mean(rank), p_same_rank = mean(rank == rank_real)),
                   by = .(tier, t, design, caller)], by = c("tier", "t", "caller"))
setorder(cbs, tier, t, design, rank_real)
fwrite(cbs, file.path(OUT, "caller_best.tsv"), sep = "\t")

# real sweep: rank stability of config F1 / precision across thresholds
sw = rbindlist(lapply(TIERS, function(tt) {
  x = dcast(real[tier == tt], config ~ t, value.var = c("F_real", "P_real"))
  pr = combn(colnames(pos), 2, simplify = FALSE)
  rbindlist(lapply(pr, function(q) data.table(tier = tt, t_a = q[1], t_b = q[2],
    rho_f1 = cor(x[[paste0("F_real_", q[1])]], x[[paste0("F_real_", q[2])]], method = "spearman", use = "complete.obs"),
    rho_precision = cor(x[[paste0("P_real_", q[1])]], x[[paste0("P_real_", q[2])]], method = "spearman", use = "complete.obs"))))
}))
fwrite(sw, file.path(OUT, "sweep_rank_stability.tsv"), sep = "\t")

print(sw)
print(rk[, .(rho_f1 = round(median(rho_f1), 3), rho_prec = round(median(rho_precision, na.rm = TRUE), 3),
             rho_prec_lo = round(quantile(rho_precision, 0.025, na.rm = TRUE), 3)), by = .(tier, t, design)])
print(cbs[design == "resolved_matched", .(tier, t, caller, config_real, F_real = round(F_real, 3), P_real = round(P_real, 3),
                                          P_rand = round(P_rand_mean, 3), rank_real, rank_rand = round(rank_rand_mean, 2),
                                          p_same_rank = round(p_same_rank, 2))])
print(desc[, .(n_calls = round(median(n_calls)), frac_not_in_table = round(median(frac_not_in_table), 3),
               frac_weak_ntc = round(median(frac_weak_ntc), 3)), by = .(tier, caller)])

# ---- figure ----------------------------------------------------------------------------
TL = c(arm = "Arm", sample = "Sample", celltype = "Cell type")
tlab = setNames(sprintf("Positive: >= %d%% depletion", 100 * TS), colnames(pos))
cols = cmap_method[std_caller(unique(long$caller))]
ro[, `:=`(Caller = std_caller(caller), Tier = factor(TL[tier], TL), T = as.numeric(sub("t", "", t)))]
pa = ggplot(melt(ro, id.vars = c("Tier", "Caller", "T"), measure.vars = c("F_real", "P_real", "R_real")),
            aes(T, value, colour = Caller)) + geom_line() + geom_point(size = 1.8) +
  facet_grid(factor(variable, c("F_real", "P_real", "R_real"), c("F1", "Precision", "Recall")) ~ Tier, scales = "free_y") +
  scale_colour_manual(values = cols) + scale_x_continuous(breaks = 100 * TS) +
  labs(x = "Depletion threshold for positives (%)", y = "Best reported configuration per caller", colour = NULL) +
  theme_bw(11) + theme(legend.position = "bottom")
cs[, `:=`(Caller = std_caller(caller), Tier = factor(TL[tier], TL), Thr = factor(tlab[t], tlab))]
pb = ggplot(cs[design == "resolved_matched"], aes(P_mean, P_real, colour = Caller)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = P_lo, xmax = P_hi), height = 0, alpha = 0.35) + geom_point(size = 1.5) +
  facet_grid(Thr ~ Tier, scales = "free") + scale_colour_manual(values = cols) +
  labs(x = "Precision, labels shuffled within resolved sites (matched strata)", y = "Precision, real labels") +
  theme_bw(10) + theme(legend.position = "none")
rk[, `:=`(Tier = factor(TL[tier], TL), T = factor(as.numeric(sub("t", "", t))),
          Design = factor(design, c("resolved_matched", "resolved_permuted", "universe_matched"),
                          c("Shuffled, matched", "Shuffled, unmatched", "Random sites, matched")))]
pc = ggplot(melt(rk, id.vars = c("Tier", "T", "Design"), measure.vars = c("rho_precision", "rho_f1")),
            aes(T, value, fill = Design)) + geom_boxplot(outlier.size = 0.3) +
  facet_grid(factor(variable, c("rho_precision", "rho_f1"), c("Precision", "F1")) ~ Tier) +
  geom_hline(yintercept = 0, colour = "grey50") +
  labs(x = "Depletion threshold (%)", y = "Spearman rho across configurations\n(real vs one random draw)", fill = NULL) +
  theme_bw(10) + theme(legend.position = "bottom")
p = plot_grid(pa, plot_grid(pb, pc, nrow = 1, labels = c("b", "c")), ncol = 1, labels = c("a", ""), rel_heights = c(1, 1.1))
ggsave(file.path(OUT, "sweep_infneg.pdf"), p, width = 16, height = 16)
ggsave(file.path(OUT, "sweep_infneg.png"), p, width = 16, height = 16, dpi = 85)
msg("done -> ", OUT)
