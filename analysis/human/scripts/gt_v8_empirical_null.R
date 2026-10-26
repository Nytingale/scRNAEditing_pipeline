#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_empirical_null.R   (2026-10-17)
#
# Empirical (decoy) null for the v8 KO positives: the positive rule is applied unchanged to the 10
# mismatch types ADAR cannot make (not A>G, not T>C; RNA read orientation) at v8-eligible positions,
# and their pass rate estimates how many A>G positives would pass by chance. Universe, strata and
# estimator: gt_v8_decoy_common.R.
#   rules           current   a covered KO arm (gate rule) with KO edit <= (1 - t) x NTC edit (rule 4)
#                   binom     current + pbinom(KO alt, KO cov, NTC edit) <= 0.05 in that arm
#                   adar1     current, ADAR1 KO arm only
# Parity: A>G positives under `current` must equal gt_eval_v8/truth/gt_v8_truth_counts.tsv (asserted).
#
# OUTPUT: gt_eval_v8/empirical_null/{empirical_null_fdr,empirical_null_strata}.tsv
# =============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "gt_v8_decoy_common.R"))
dir.create(DECOY_OUT, showWarnings = FALSE)
ALPHA = 0.05
v8n = fread(file.path(V8, "truth/gt_v8_truth_counts.tsv"))

res = list(); strata = list()
for (g in names(GATES)) {
  d = load_v8_universe(g)
  stopifnot(d[grp == "AG", .N] == v8n[gate == g & t == "t25", eligible])
  covered = lapply(ARMS, function(a) {
    cov = d[[paste0(a, "_cov")]]
    x = d$ntc_cov >= 10 & if (g == "expected") cov * d$ntc_edit >= 3 else cov >= 10
    x[is.na(x)] = FALSE; x })
  pass = list()
  for (tn in names(DS)) {
    pt = lapply(names(ARMS), function(k) {
      a = ARMS[[k]]; ke = d[[paste0(a, "_ke")]]
      x = covered[[k]] & !is.na(ke) & ke <= (1 - DS[[tn]]) * d$ntc_edit; x[is.na(x)] = FALSE; x })
    pb = lapply(seq_along(ARMS), function(i) pt[[i]] & d[[paste0(ARMS[[i]], "_p")]] <= ALPHA)
    pass[[paste(tn, "current")]] = pt[[1]] | pt[[2]]
    pass[[paste(tn, "binom")]] = pb[[1]] | pb[[2]]
    pass[[paste(tn, "adar1")]] = pt[[1]]
    stopifnot(sum(pass[[paste(tn, "current")]][d$grp == "AG"]) == v8n[gate == g & t == tn, positive])
  }
  r = decoy_fdr(d, do.call(cbind, pass))
  r[, c("t", "rule") := tstrsplit(rule, " ", fixed = TRUE)]
  res[[g]] = r[, gate := g]
  strata[[g]] = strata_table(d)[, gate := g]
  message(g, ": ", d[grp == "AG", .N], " A>G, ", d[grp == "TC", .N], " T>C, ", d[grp == "decoy", .N], " decoy sites")
}
res = rbindlist(res)
setcolorder(res, c("gate", "t", "rule"))
fwrite(res, file.path(DECOY_OUT, "empirical_null_fdr.tsv"), sep = "\t")
fwrite(rbindlist(strata), file.path(DECOY_OUT, "empirical_null_strata.tsv"), sep = "\t")
print(res[, .(gate, t, rule, ag_positives, ag = round(ag_rate, 3), tc = round(tc_rate, 3), decoy = round(decoy_rate, 3),
              expected = round(expected_chance), fdr = round(fdr, 2), ci = sprintf("%.2f-%.2f", fdr_lo, fdr_hi),
              above = round(above_chance))])
