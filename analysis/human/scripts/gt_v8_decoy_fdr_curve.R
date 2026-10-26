#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_decoy_fdr_curve.R   (2026-10-17)
#
# Decoy-calibrated positives without the m (expected KO reads) gate. A site is positive at (p, t) when
# a KO arm has KO edit <= (1 - t) x NTC edit and pbinom(KO alt, KO cov, NTC edit) <= p in that arm,
# i.e. the site is ranked by the smaller of its two arm p-values (no multiple-testing correction; the
# decoys set the threshold). Any arm with KO cov > 0 counts. p = 1 is the effect size alone.
# Universe, strata and estimator: gt_v8_decoy_common.R (sites of each cover gate's v8 universe).
#   overdispersion  per arm, share of decoy arms (KO cov > 0) with p <= nominal p; with a binomial
#                   null and discrete counts this share is at most the nominal p
#
# OUTPUT: gt_eval_v8/empirical_null/{decoy_fdr_curve,decoy_overdispersion}.tsv
# =============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "gt_v8_decoy_common.R"))
dir.create(DECOY_OUT, showWarnings = FALSE)
P_THR = c(1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 0.01, 0.025, 0.05, 0.1, 0.2, 1)

curve = list(); od = list()
for (g in names(GATES)) {
  d = load_v8_universe(g)
  pass = list()
  for (tn in names(DS)) for (p in P_THR) {
    arm_ok = lapply(ARMS, function(a) {
      ke = d[[paste0(a, "_ke")]]; pv = d[[paste0(a, "_p")]]
      x = !is.na(ke) & ke <= (1 - DS[[tn]]) * d$ntc_edit & pv <= p; x[is.na(x)] = FALSE; x })
    pass[[paste(tn, p)]] = arm_ok[[1]] | arm_ok[[2]]
  }
  r = decoy_fdr(d, do.call(cbind, pass))
  r[, c("t", "p_thr") := tstrsplit(rule, " ", fixed = TRUE)][, p_thr := as.numeric(p_thr)]
  curve[[g]] = r[, `:=`(gate = g, rule = NULL)]
  od[[g]] = rbindlist(lapply(names(ARMS), function(k) {
    a = ARMS[[k]]; x = d[get(paste0(a, "_cov")) > 0]
    rbindlist(lapply(P_THR[P_THR < 1], function(p)
      x[, .(n_arms = .N, share = mean(get(paste0(a, "_p")) <= p)), by = grp][, `:=`(arm = k, nominal_p = p)]))
  }))[, gate := g]
}
curve = rbindlist(curve)
setcolorder(curve, c("gate", "t", "p_thr"))
fwrite(curve, file.path(DECOY_OUT, "decoy_fdr_curve.tsv"), sep = "\t")
od = rbindlist(od)
fwrite(od, file.path(DECOY_OUT, "decoy_overdispersion.tsv"), sep = "\t")
print(curve[, .(gate, t, p_thr, ag_positives, expected = round(expected_chance), fdr = round(fdr, 2),
                ci = sprintf("%.2f-%.2f", fdr_lo, fdr_hi), above = round(above_chance))], nrows = 200)
print(dcast(od[gate == "expected"], arm + nominal_p ~ grp, value.var = "share"), digits = 3)
