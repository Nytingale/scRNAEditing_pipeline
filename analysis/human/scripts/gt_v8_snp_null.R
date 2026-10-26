#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_snp_null.R   (2026-10-17)
#
# SNP-calibrated null for A>G positives. PGP1 SNVs cannot respond to an ADAR KO (all arms are PGP1
# cells), so their NTC vs KO allele-fraction differences are sampling plus technical noise.
#   SNP set     12-type gate sites of each v8 universe at PGP1 SNV positions (removed from the truth by the
#               PGP1 filter), nearby-mismatch pass (computed before PGP1 removal), positions with one type
#   noise model KO alt ~ beta-binomial(KO cov, mean = NTC allele fraction, overdispersion rho), rho fitted
#               per KO arm by maximum likelihood on heterozygous SNPs (0.1 <= NTC AF < 0.9, KO cov > 0);
#               95% CI from 200 bootstrap draws of the SNPs. Models: binomial (rho = 0, the current test),
#               beta_binomial, bb_ntc (main: rho + 1 / NTC cov, adding the sampling error of the NTC level)
#   positive    an arm with KO edit <= (1 - t) x NTC edit and beta-binomial lower-tail p <= p_thr
#               (rank by the smaller arm p, as gt_v8_decoy_fdr_curve.R)
#   expected chance positives = sum over A>G sites of P(site passes | null model), arms independent;
#               FDR upper bound = expected / observed (assumes every A>G site is unaffected by the KO);
#               range from rho at its CI limits
#   checks      the same expectation vs observed passes for decoys (10 non-A>G non-T>C types, by NTC AF
#               bin), held-out SNPs (rho fitted on odd chromosomes, tested on even, and vice versa) and
#               T>C (positive control: should exceed); rho of A>G / T>C SNPs vs other SNPs (likelihood ratio)
#
# OUTPUT: gt_eval_v8/snp_null/{snp_rho,snp_rho_by_cov,snp_null_calibration,snp_null_fdr}.tsv,
#         snp_null_ag_sites.tsv.gz (per A>G site and arm: beta-binomial p)
# =============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "gt_v8_decoy_common.R"))
OUT = file.path(V8, "snp_null"); dir.create(OUT, showWarnings = FALSE)
P_THR = c(1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 0.01, 0.025, 0.05, 0.1, 0.2, 1)
HET = c(0.1, 0.9)
N_BOOT = 200
RHO_MIN = 1e-8
MODELS = list(binomial = list(fit = FALSE, ntc = FALSE), beta_binomial = list(fit = TRUE, ntc = FALSE),
              bb_ntc = list(fit = TRUE, ntc = TRUE))
MAIN = "bb_ntc"

load_snp_universe = function(g) {
  flags = fread(file.path(V8, "truth/support/truth12_site_flags.tsv.gz"),
                select = c("gate", "chr", "pos", "strand", "type", "pgp1", "nomm_pass"))[gate == g]
  d = fread(file.path(PM, "gt_eval_v6/spectrum12", GATES[[g]], "truth12_sites.tsv.gz"))
  d = merge(d, flags, by = c("chr", "pos", "strand", "type"))[pgp1 == TRUE & nomm_pass == TRUE]
  d[, site := paste0(chr, ":", pos)]
  d = d[!site %chin% d[, .N, by = site][N > 1, site]]
  d[, grp := "SNP"][, snp_class := fifelse(type %in% c("AG", "TC"), "AG_TC", "other")]
  annotate_arms(d)
}

lbb = function(k, n, mu, rho) {
  s = 1 / rho - 1; a = mu * s; b = (1 - mu) * s
  lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
}
lbinom = function(k, n, mu) dbinom(k, n, mu, log = TRUE)
# per-site overdispersion; ntc = TRUE adds the sampling error of the NTC level (Var(NTC AF) / (mu (1 - mu))
# = 1 / NTC cov), so KO alt ~ beta-binomial(KO cov, NTC AF, rho + 1 / NTC cov) by moment matching
rho_of = function(r, ntc_cov, ntc) pmin(r + if (ntc) 1 / ntc_cov else 0 * ntc_cov, 0.5)
fit_rho = function(alt, cov, mu, ntc_cov, ntc) {
  mu = pmin(pmax(mu, 1e-4), 1 - 1e-4)
  optimize(function(r) -sum(lbb(alt, cov, mu, rho_of(r, ntc_cov, ntc))), c(RHO_MIN, 0.5))
}

# per arm: observed pass (logical) and null pass probability, sites x (t, p_thr); rho per site
arm_eval = function(alt, cov, mu, rho) {
  combos = CJ(t = names(DS), p_thr = P_THR, sorted = FALSE)
  n = length(cov); obs = matrix(FALSE, n, nrow(combos)); pnull = matrix(0, n, nrow(combos))
  mu = pmin(pmax(mu, 1e-4), 1 - 1e-4)
  for (i in seq_len(n)) {
    if (is.na(cov[i]) || cov[i] <= 0) next
    ks = 0:min(cov[i], floor((1 - min(DS)) * mu[i] * cov[i]) + 1)
    lp = if (rho[i] <= RHO_MIN) lbinom(ks, cov[i], mu[i]) else lbb(ks, cov[i], mu[i], rho[i])
    cdf = cumsum(exp(lp))
    ke = rnd4(ks / cov[i])
    for (j in seq_len(nrow(combos))) {
      okk = ke <= (1 - DS[[combos$t[j]]]) * rnd4(mu[i]) & cdf <= combos$p_thr[j]
      if (any(okk)) pnull[i, j] = cdf[max(which(okk))]
      if (alt[i] <= max(ks)) obs[i, j] = okk[alt[i] + 1]
    }
  }
  list(obs = obs, pnull = pnull, combos = combos)
}

# site level over both arms: observed positives and expected chance positives per (t, p_thr)
site_eval = function(d, rho, ntc) {
  r = lapply(names(ARMS), function(k) {
    a = ARMS[[k]]
    arm_eval(d[[paste0(a, "_alt")]], d[[paste0(a, "_cov")]], d$ntc_edit, rho_of(rho[[k]], d$ntc_cov, ntc))
  })
  list(obs = r[[1]]$obs | r[[2]]$obs, pnull = 1 - (1 - r[[1]]$pnull) * (1 - r[[2]]$pnull), combos = r[[1]]$combos)
}

summarise = function(se, idx, label) {
  data.table(se$combos, group = label, n_sites = sum(idx),
             observed = colSums(se$obs[idx, , drop = FALSE]), expected = colSums(se$pnull[idx, , drop = FALSE]))
}

rho_tab = list(); rho_cov = list(); calib = list(); fdr = list(); sites = list()
for (g in names(GATES)) {
  d = load_v8_universe(g)
  s = load_snp_universe(g)
  s[, odd := suppressWarnings(as.integer(sub("chr", "", chr))) %% 2 == 1][is.na(odd), odd := FALSE]
  het = s[ntc_edit >= HET[1] & ntc_edit < HET[2]]
  message(g, ": SNPs ", nrow(s), " (heterozygous ", nrow(het), "); A>G ", sum(d$grp == "AG"),
          ", T>C ", sum(d$grp == "TC"), ", decoys ", sum(d$grp == "decoy"))


  # ---- overdispersion per arm and model ---------------------------------------------------------
  rho = list(); rho_ci = list(); rho_cv = list()
  for (m in names(MODELS)) {
    ntc = MODELS[[m]]$ntc
    rho[[m]] = list(); rho_ci[[m]] = list(); rho_cv[[m]] = list()
    for (k in names(ARMS)) {
      if (!MODELS[[m]]$fit) {
        rho[[m]][[k]] = 0; rho_ci[[m]][[k]] = c(0, 0); rho_cv[[m]][[k]] = c(odd = 0, even = 0); next
      }
      set.seed(20261017)
      a = ARMS[[k]]; h = het[get(paste0(a, "_cov")) > 0]
      alt = h[[paste0(a, "_alt")]]; cov = h[[paste0(a, "_cov")]]; mu = h$ntc_edit; nc = h$ntc_cov
      fr = function(j) fit_rho(alt[j], cov[j], mu[j], nc[j], ntc)
      f = fr(seq_along(alt)); rho[[m]][[k]] = f$minimum
      rho_ci[[m]][[k]] = quantile(vapply(seq_len(N_BOOT), function(i) fr(sample.int(length(alt), replace = TRUE))$minimum, 0),
                                  c(0.025, 0.975))
      llb = sum(lbb(alt, cov, pmin(pmax(mu, 1e-4), 1 - 1e-4), rho_of(RHO_MIN, nc, ntc)))
      fc = lapply(c("AG_TC", "other"), function(cl) fr(which(h$snp_class == cl)))
      rho_cv[[m]][[k]] = c(odd = fr(which(h$odd))$minimum, even = fr(which(!h$odd))$minimum)
      rho_tab[[paste(g, m, k)]] = data.table(gate = g, model = m, arm = k, n_het_snps = nrow(h),
        median_ko_cov = median(cov), median_ntc_cov = median(nc),
        rho = f$minimum, rho_lo = rho_ci[[m]][[k]][1], rho_hi = rho_ci[[m]][[k]][2],
        lr_vs_rho0 = 2 * (-f$objective - llb), rho_ag_tc = fc[[1]]$minimum, rho_other = fc[[2]]$minimum,
        n_ag_tc = sum(h$snp_class == "AG_TC"),
        p_ag_tc_vs_other = pchisq(2 * (f$objective - fc[[1]]$objective - fc[[2]]$objective), 1, lower.tail = FALSE),
        rho_odd_chr = rho_cv[[m]][[k]][["odd"]], rho_even_chr = rho_cv[[m]][[k]][["even"]])
      h[, cb := cut(cov, c(0, 5, 10, 20, 50, Inf), c("1-4", "5-9", "10-19", "20-49", ">=50"), right = FALSE)]
      rho_cov[[paste(g, m, k)]] = h[, .(n = .N, rho = fit_rho(get(paste0(a, "_alt")), get(paste0(a, "_cov")),
                                                              ntc_edit, ntc_cov, ntc)$minimum),
                                    by = cb][, `:=`(gate = g, model = m, arm = k)]
    }
  }

  # ---- expected vs observed ---------------------------------------------------------------------
  ag = d$grp == "AG"
  for (m in names(MODELS)) {
    ntc = MODELS[[m]]$ntc
    se = site_eval(d, rho[[m]], ntc)
    cl = list(summarise(se, d$grp == "AG", "AG"), summarise(se, d$grp == "TC", "TC"),
              summarise(se, d$grp == "decoy", "decoy"))
    cl = lapply(cl, function(x) x[, s_af := "all"])
    for (af in AF_LAB) if (any(d$grp == "decoy" & d$s_af == af))
      cl[[length(cl) + 1]] = summarise(se, d$grp == "decoy" & d$s_af == af, "decoy")[, s_af := af]
    # held-out heterozygous SNPs: rho from the other chromosome half
    ho = rbindlist(lapply(c(TRUE, FALSE), function(o) {
      r = lapply(rho_cv[[m]], function(v) v[[if (o) "even" else "odd"]])
      summarise(site_eval(het[odd == o], r, ntc), rep(TRUE, sum(het$odd == o)), "SNP_heldout")
    }))
    cl[[length(cl) + 1]] = ho[, .(n_sites = sum(n_sites), observed = sum(observed), expected = sum(expected)),
                              by = .(t, p_thr, group)][, s_af := "all"]
    calib[[paste(g, m)]] = rbindlist(cl, use.names = TRUE)[, `:=`(gate = g, model = m)]

    lo = colSums(site_eval(d[ag], lapply(rho_ci[[m]], `[`, 1), ntc)$pnull)
    hi = colSums(site_eval(d[ag], lapply(rho_ci[[m]], `[`, 2), ntc)$pnull)
    f = data.table(se$combos, gate = g, model = m, ag_sites = sum(ag), ag_positives = colSums(se$obs[ag, ]),
                   expected = colSums(se$pnull[ag, ]), expected_rho_lo = lo, expected_rho_hi = hi)
    f[, `:=`(fdr = expected / pmax(ag_positives, 1), fdr_lo = pmin(lo, hi) / pmax(ag_positives, 1),
             fdr_hi = pmax(lo, hi) / pmax(ag_positives, 1))]
    fdr[[paste(g, m)]] = f
  }

  # per A>G site and arm: lower-tail p under the main model
  x = d[ag, c("site", "strand", "ntc_cov", "ntc_alt", "ntc_edit", paste0(rep(ARMS, each = 2), c("_cov", "_alt"))), with = FALSE]
  for (k in names(ARMS)) {
    a = ARMS[[k]]; cv = x[[paste0(a, "_cov")]]; al = x[[paste0(a, "_alt")]]
    mu = pmin(pmax(x$ntc_edit, 1e-4), 1 - 1e-4); rr = rho_of(rho[[MAIN]][[k]], x$ntc_cov, MODELS[[MAIN]]$ntc)
    p = vapply(seq_len(nrow(x)), function(i) if (is.na(cv[i]) || cv[i] <= 0) NA_real_ else
      sum(exp(lbb(0:al[i], cv[i], mu[i], rr[i]))), 0)
    set(x, j = paste0(a, "_p"), value = pmin(p, 1))
  }
  sites[[g]] = x[, `:=`(gate = g, model = MAIN)]
}

rho_tab = rbindlist(rho_tab); fwrite(rho_tab, file.path(OUT, "snp_rho.tsv"), sep = "\t")
rho_cov = rbindlist(rho_cov); fwrite(rho_cov, file.path(OUT, "snp_rho_by_cov.tsv"), sep = "\t")
calib = rbindlist(calib, use.names = TRUE)
setcolorder(calib, c("gate", "model", "group", "s_af", "t", "p_thr"))
fwrite(calib, file.path(OUT, "snp_null_calibration.tsv"), sep = "\t")
fdr = rbindlist(fdr); setcolorder(fdr, c("gate", "model", "t", "p_thr"))
fwrite(fdr, file.path(OUT, "snp_null_fdr.tsv"), sep = "\t")
fwrite(rbindlist(sites), file.path(OUT, "snp_null_ag_sites.tsv.gz"), sep = "\t")

print(rho_tab, digits = 3)
print(dcast(rho_cov, gate + model + arm ~ cb, value.var = "rho"), digits = 3)
print(dcast(calib[t == "t50" & s_af == "all" & group != "AG" & p_thr %in% c(1e-3, 2.5e-3, 0.01, 0.05, 1)],
            gate + group + p_thr ~ model, value.var = c("observed", "expected")), digits = 3, nrows = 200)
print(dcast(calib[gate == "expected" & group == "decoy" & s_af != "all" & t == "t50" & p_thr %in% c(0.01, 0.05)],
            s_af + p_thr ~ model, value.var = c("observed", "expected")), digits = 3)
print(fdr[model == MAIN, .(gate, t, p_thr, ag_positives, expected = round(expected, 1), fdr = round(fdr, 2),
                           range = sprintf("%.2f-%.2f", fdr_lo, fdr_hi))], nrows = 200)
