#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_m_null_check_plot.R   (2026-10-18)
#
# Review response, Set A of HANDOFF_2026-10-17 4c-DECISION (m = 3 kept as an eligibility gate). Source data
# for Figure 2 e-f (GTDEF_gtdefinition_plot.R, GTDEF_v8_informative_ko) and a standalone copy of those panels.
#   e  share of KO arms with 0 edited (alt) reads, binned by m = KO coverage x NTC allele fraction (bins
#      [k - 0.5, k + 0.5), plotted at k = 1..8), against e^-m: non-AG/TC sites (cannot respond to an ADAR KO) = the 10
#      mismatch types that are neither A>G nor T>C (RNA orientation); and A-to-G candidate sites (RNA A>G, both
#      genomic strands; chance zeros + KO effect). PGP1 SNV positions are removed from both. One row per site and KO (ADAR1 KO, ADAR2 KO) with KO
#      coverage >= 1).
#      Plotted (user 2026-10-18): binomial GAM of zero ~ s(m) per group and KO on the sites (0 < m <= 8.5;
#      each KO separately, so every site is counted once per fit),
#      fit and 95% band -> GTDEF_v8_m_null_fit.tsv. The binned table (GTDEF_v8_m_null_zero.tsv) is kept for
#      the counts quoted in the text.
#   (not plotted; user 2026-10-18: quote only if asked) share of A-to-G zero-read KO arms beyond chance
#      = 1 - non-AG/TC rate / A-to-G rate per m bin; 95% interval from 2,000 binomial draws of both rates; bins with
#      >= 10 A-to-G zero-read arms -> GTDEF_v8_m_null_excess.tsv
# Site universe: the v8 rules with the gate relaxed to NTC coverage >= 10 and KO coverage >= 1
# (gate depth1; run_gt_v8_fisher.sh), i.e. every site that passes the NTC rules and has a KO read.
#
# INPUT : human/pm/gt_eval_v6/spectrum12/eval_depth1/truth12_sites.tsv.gz,
#         human/pm/gt_eval_v8/fisher/truth/{truth12_site_flags.tsv.gz, v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz}
# OUTPUT: report/pm_gt/GTDEF_v8_m_null_zero.tsv, report/pm_gt/GTDEF_v8_m_null_excess.tsv,
#         figures/pm_gt/GTDEF_v8_m_null_check.svg
# Usage : Rscript GTDEF_v8_m_null_check_plot.R
# ============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "GTDEF_common.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
gtdef_set_dataset("pm_gt")
SCRIPT = caller_script()
PM = file.path(REPO_DIR, "human", "pm")
FT = file.path(PM, "gt_eval_v8", "fisher")
M_K = 1:8
M_BR = c(M_K - 0.5, max(M_K) + 0.5)
N_DRAW = 2000; MIN_AG_ZEROS = 10
G_DEC = "Non-AG/TC sites (10 mismatch types)"
G_AG = "A-to-G candidate sites"

t12 = fread(file.path(PM, "gt_eval_v6/spectrum12/eval_depth1/truth12_sites.tsv.gz"))
flags = fread(file.path(FT, "truth/truth12_site_flags.tsv.gz"),
              select = c("gate", "chr", "pos", "strand", "type", "pgp1", "nomm_pass", "nomm_pass_pgp1"))[gate == "depth1"]
truth = fread(file.path(FT, "truth/v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz"), select = c("gate", "site"))[gate == "depth1"]
t12 = merge(t12, flags[, gate := NULL], by = c("chr", "pos", "strand", "type"))
t12[, site := paste0(chr, ":", pos)]
t12[, ntc_af := ntc_alt / ntc_cov]
one_type = function(x) x[!site %chin% x[, .N, by = site][N > 1, site]]

u = t12[pgp1 == FALSE & nomm_pass_pgp1 == TRUE]
dec = one_type(u[!type %in% c("AG", "TC")])[, grp := G_DEC]
ag = u[type == "AG" & site %chin% truth$site][, grp := G_AG]
stopifnot(nrow(ag) == nrow(truth))

arms_long = function(x) rbindlist(lapply(c(ADAR1 = "adar1", ADAR2 = "adarb1"), function(a)
  x[, .(grp, site, ntc_af, cov = get(paste0(a, "_cov")), alt = get(paste0(a, "_alt")))]), idcol = "arm")[cov >= 1]
nl = arms_long(rbind(dec, ag))
nl[, `:=`(m = cov * ntc_af, zero = alt == 0, p_binom = (1 - ntc_af)^cov)]
fit_grid = seq(0.1, max(M_K), 0.05)
zf = rbindlist(lapply(split(nl[m > 0 & m <= max(M_K) + 0.5], by = c("grp", "arm")), function(x) {
  g = mgcv::gam(zero ~ s(m, k = 10), family = binomial, data = x, method = "REML")
  p = predict(g, data.frame(m = fit_grid), se.fit = TRUE)
  data.table(grp = x$grp[1], arm = x$arm[1], m = fit_grid, fit = plogis(p$fit), lo = plogis(p$fit - 1.96 * p$se.fit),
             hi = plogis(p$fit + 1.96 * p$se.fit), n_sites = nrow(x))
}))
fwrite(zf, file.path(REPORT_DIR, "GTDEF_v8_m_null_fit.tsv"), sep = "\t")
nl[, mbin := cut(m, M_BR, right = FALSE)]
nl = nl[!is.na(mbin)]
nl[, m_k := M_K[as.integer(mbin)]]
z = nl[, .(n = .N, zeros = sum(zero), m_mean = mean(m), p_binom = mean(p_binom), p_poisson = mean(exp(-m))),
       by = .(grp, arm, mbin, m_mid = m_k)]
z[, `:=`(obs = zeros / n, lo = qbeta(0.025, zeros, n - zeros + 1), hi = qbeta(0.975, zeros + 1, n - zeros))]
z[zeros == 0, lo := 0]
z = z[n >= 20]
setorder(z, grp, arm, mbin)
fwrite(z, file.path(REPORT_DIR, "GTDEF_v8_m_null_zero.tsv"), sep = "\t")

w = merge(z[grp == G_AG, .(arm, mbin, m_mid, n_ag = n, z_ag = zeros)], z[grp == G_DEC, .(arm, mbin, n_null = n, z_null = zeros)],
          by = c("arm", "mbin"))
w = w[z_ag >= MIN_AG_ZEROS]
set.seed(20261018)
w[, c("excess", "lo", "hi") := {
  e = 1 - (z_null / n_null) / (z_ag / n_ag)
  d = vapply(seq_len(.N), function(i) {
    ra = rbinom(N_DRAW, n_ag[i], z_ag[i] / n_ag[i]) / n_ag[i]; rn = rbinom(N_DRAW, n_null[i], z_null[i] / n_null[i]) / n_null[i]
    q = quantile(1 - rn / pmax(ra, 1e-9), c(0.025, 0.975)); c(q[[1]], q[[2]]) }, numeric(2))
  list(e, pmax(d[1, ], 0), pmin(d[2, ], 1)) }]
fwrite(w, file.path(REPORT_DIR, "GTDEF_v8_m_null_excess.tsv"), sep = "\t")
print(z[, .(grp, mbin, n, zeros, obs = round(obs, 3), e_m = round(p_poisson, 3))], nrows = 200)
print(w, digits = 3)

# standalone copy of Figure 2 e-f (the same builders as GTDEF_gtdefinition_plot.R)
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])),
                 "GTDEF_v8_m_null_panels.R"))
fig = m_null_zero_panel(zf, 13) + theme(legend.position = "right")
save_figure(fig, "GTDEF_v8_m_null_check", 15, 5.6, file.path(REPORT_DIR, "GTDEF_v8_m_null_fit.tsv"),
            data = zf, script = SCRIPT)
