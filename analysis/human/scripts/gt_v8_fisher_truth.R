#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_fisher_truth.R   (2026-10-18)
#
# Set B of HANDOFF_2026-10-17 4c-DECISION: KO positives by a one-sided Fisher's exact test instead of the
# m = 3 expected-read gate. Universe = the v8 truth built with gate depth1 (NTC cov >= 10, KO cov >= 1 per
# arm; all other v8 rules and cleaning unchanged; run_gt_v8_fisher.sh stages eval + truth).
#   Fisher p per arm  one-sided (KO editing lower than merged NTC): phyper(KO alt, KO alt + NTC alt,
#                     KO ref + NTC ref, KO cov) = fisher.test(alternative = "less") on the 2 x 2 read table
#   depth1            effect size only: an arm with KO cov >= 1 and KO edit <= (1 - t) x NTC edit (parity with
#                     the depth1 truth asserted)
#   fisher_bonf       depth1 rule in an arm whose Fisher p <= 0.025 (Bonferroni for the two arms)
#   fisher_bh         per-site p = min(1, 2 x smaller arm p); BH across the eligible A>G sites of the universe;
#                     an arm passes if 2 x its p <= the largest per-site p with q <= 0.10 and it meets the
#                     depth1 rule (BH once, independent of t)
# Informative negatives are unchanged (rule 5); positives failing the test become unresolved.
# Chance positives per rule and t:
#   decoy   gt_v8_decoy_common.R decoy_fdr() on the depth1 12-type universe (10 non-A>G non-T>C types);
#           the BH rule is applied to decoys with the A>G cut-off
#   SNP     KO alt ~ beta-binomial(KO cov, NTC AF, rho + 1 / NTC cov) (gt_v8_snp_null.R model bb_ntc), rho per
#           arm fitted on the depth1 heterozygous PGP1 SNPs; expected = sum over A>G sites of P(site passes)
#
# INPUT : gt_eval_v6/spectrum12/eval_depth1/truth12_sites.tsv.gz,
#         gt_eval_v8/fisher/truth/{truth12_site_flags.tsv.gz, v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz}
# OUTPUT: gt_eval_v8/fisher/truth_fisher/gt_v7_truth_sites.tsv.gz (gates depth1, fisher_bonf, fisher_bh;
#         gt_score_v7.R GT_V7_TRUTH), fisher_truth_counts.tsv, fisher_site_p.tsv.gz, fisher_chance.tsv,
#         fisher_snp_rho.tsv
# =============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "gt_v8_decoy_common.R"))
FT = file.path(V8, "fisher")
OUT = file.path(FT, "truth_fisher"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
G0 = "depth1"
ALPHA_BONF = 0.025; Q_BH = 0.10
HET = c(0.1, 0.9); RHO_MIN = 1e-8; N_BOOT = 200

fisher_less = function(ka, kc, na, nc) fifelse(!is.na(kc) & kc > 0, phyper(ka, ka + na, (kc - ka) + (nc - na), kc), NA_real_)

# ---- universe ---------------------------------------------------------------------------------------
truth = fread(file.path(FT, "truth/v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz"))[gate == G0]
flags = fread(file.path(FT, "truth/truth12_site_flags.tsv.gz"),
              select = c("gate", "chr", "pos", "strand", "type", "pgp1", "nomm_pass", "nomm_pass_pgp1"))[gate == G0][, gate := NULL]
t12 = fread(file.path(PM, "gt_eval_v6/spectrum12/eval_depth1/truth12_sites.tsv.gz"))
t12 = merge(t12, flags, by = c("chr", "pos", "strand", "type"))
t12[, site := paste0(chr, ":", pos)]

d = t12[pgp1 == FALSE & nomm_pass_pgp1 == TRUE]
d = d[type != "AG" | site %chin% truth$site]
d = d[type == "AG" | !site %chin% d[, .N, by = site][N > 1, site]]
d[, grp := fifelse(type == "AG", "AG", fifelse(type == "TC", "TC", "decoy"))]
invisible(annotate_arms(d))
stopifnot(d[grp == "AG", .N] == nrow(truth), d[grp == "AG", uniqueN(site)] == nrow(truth))
for (a in ARMS) set(d, j = paste0(a, "_fp"),
                    value = fisher_less(d[[paste0(a, "_alt")]], d[[paste0(a, "_cov")]], d$ntc_alt, d$ntc_cov))
d[, p_site := pmin(1, 2 * pmin(fcoalesce(adar1_fp, 1), fcoalesce(adarb1_fp, 1)))]
ag = d$grp == "AG"
q = p.adjust(d$p_site[ag], "BH")
P_BH = if (any(q <= Q_BH)) max(d$p_site[ag][q <= Q_BH]) else 0
message(G0, ": A>G ", sum(ag), ", T>C ", sum(d$grp == "TC"), ", decoys ", sum(d$grp == "decoy"),
        "; BH q <= ", Q_BH, ": ", sum(q <= Q_BH), " sites, per-site p cut-off ", signif(P_BH, 3))

# arm meets the effect-size rule at t (KO cov >= 1; edits rounded as classify())
arm_drop = function(a, t) { cv = d[[paste0(a, "_cov")]]; !is.na(cv) & cv > 0 & d[[paste0(a, "_ke")]] <= (1 - t) * d$ntc_edit }
arm_p = function(a) fcoalesce(d[[paste0(a, "_fp")]], 1)
RULES = list(depth1 = function(a) rep(TRUE, nrow(d)),
             fisher_bonf = function(a) arm_p(a) <= ALPHA_BONF,
             fisher_bh = function(a) 2 * arm_p(a) <= P_BH)
pass = list()
for (r in names(RULES)) for (tn in names(DS))
  pass[[paste(r, tn, sep = "__")]] = Reduce(`|`, lapply(ARMS, function(a) arm_drop(a, DS[[tn]]) & RULES[[r]](a)))
pass = as.data.table(pass)

# ---- truth tables -------------------------------------------------------------------------------
x = merge(truth, cbind(d[ag, .(site, read_strand = strand)], pass[ag]), by = c("site", "read_strand"))
stopifnot(nrow(x) == nrow(truth))
for (tn in names(DS)) stopifnot(identical(x[[paste0("class_", tn)]] == "positive", x[[paste0("depth1__", tn)]]))
message("parity: depth1 effect-size positives reproduce the depth1 truth at t25 / t50 / t75")
COLS = names(truth)
out = rbindlist(lapply(names(RULES), function(r) {
  y = copy(x)[, gate := r]
  for (tn in names(DS)) {
    cl = paste0("class_", tn)
    set(y, j = cl, value = fifelse(y[[paste(r, tn, sep = "__")]], "positive",
                                   fifelse(y[[cl]] == "informative_negative", "informative_negative", "unresolved")))
  }
  y[, ..COLS]
}))
setorder(out, gate, chr, pos)
fwrite(out, file.path(OUT, "gt_v7_truth_sites.tsv.gz"), sep = "\t")
cnt = rbindlist(lapply(names(DS), function(tn) { cl = paste0("class_", tn)
  out[, .(t = tn, eligible = .N, positive = sum(get(cl) == "positive"),
          informative_negative = sum(get(cl) == "informative_negative")), by = gate] }))
cnt[, bh_p_cutoff := fifelse(gate == "fisher_bh", P_BH, NA_real_)]
fwrite(cnt, file.path(OUT, "fisher_truth_counts.tsv"), sep = "\t")
print(dcast(cnt, gate ~ t, value.var = "positive"))
fwrite(d[ag, .(site, strand, ntc_cov, ntc_alt, ntc_edit, adar1_cov, adar1_alt, adar1_fp, adarb1_cov, adarb1_alt,
               adarb1_fp, p_site, q_bh = q)], file.path(OUT, "fisher_site_p.tsv.gz"), sep = "\t")

# ---- chance positives: decoys ---------------------------------------------------------------------
dec = decoy_fdr(d, pass)
dec[, c("rule", "t") := tstrsplit(rule, "__", fixed = TRUE)][, null := "decoy"]

# ---- chance positives: SNP-calibrated beta-binomial -------------------------------------------------
lbb = function(k, n, mu, rho) { s = 1 / rho - 1; lchoose(n, k) + lbeta(k + mu * s, n - k + (1 - mu) * s) - lbeta(mu * s, (1 - mu) * s) }
rho_of = function(r, ntc_cov) pmin(r + 1 / ntc_cov, 0.5)
fit_rho = function(alt, cov, mu, nc) {
  mu = pmin(pmax(mu, 1e-4), 1 - 1e-4)
  optimize(function(r) -sum(lbb(alt, cov, mu, rho_of(r, nc))), c(RHO_MIN, 0.5))$minimum
}
s = t12[pgp1 == TRUE & nomm_pass == TRUE]
s = s[!site %chin% s[, .N, by = site][N > 1, site]]
invisible(annotate_arms(s))
het = s[ntc_edit >= HET[1] & ntc_edit < HET[2]]
rho = list(); rho_tab = list()
for (k in names(ARMS)) {
  a = ARMS[[k]]; h = het[get(paste0(a, "_cov")) > 0]
  al = h[[paste0(a, "_alt")]]; cv = h[[paste0(a, "_cov")]]
  rho[[k]] = fit_rho(al, cv, h$ntc_edit, h$ntc_cov)
  set.seed(20261018)
  ci = quantile(vapply(seq_len(N_BOOT), function(i) { j = sample.int(length(al), replace = TRUE)
    fit_rho(al[j], cv[j], h$ntc_edit[j], h$ntc_cov[j]) }, 0), c(0.025, 0.975))
  rho_tab[[k]] = data.table(gate = G0, arm = k, n_het_snps = nrow(h), median_ko_cov = median(cv), rho = rho[[k]],
                            rho_lo = ci[[1]], rho_hi = ci[[2]])
}
rho_tab = rbindlist(rho_tab)
fwrite(rho_tab, file.path(OUT, "fisher_snp_rho.tsv"), sep = "\t")
print(rho_tab, digits = 3)

# P(arm passes | null) per A>G site: the pass set in KO alt is {0..kmax} (Fisher p and KO edit rise with alt)
arm_pnull = function(a, rho_a, t, rule) {
  cv = d[[paste0(a, "_cov")]][ag]; na = d$ntc_alt[ag]; nc = d$ntc_cov[ag]; mu0 = d$ntc_edit[ag]
  mu = pmin(pmax(mu0, 1e-4), 1 - 1e-4); rr = rho_of(rho_a, nc)
  vapply(seq_along(cv), function(i) {
    if (is.na(cv[i]) || cv[i] <= 0) return(0)
    ks = 0:min(cv[i], floor((1 - t) * mu0[i] * cv[i]) + 1)
    ok = rnd4(ks / cv[i]) <= (1 - t) * mu0[i]
    if (rule != "depth1") {
      fp = phyper(ks, ks + na[i], (cv[i] - ks) + (nc[i] - na[i]), cv[i])
      ok = ok & (if (rule == "fisher_bonf") fp <= ALPHA_BONF else 2 * fp <= P_BH)
    }
    if (!any(ok)) return(0)
    sum(exp(lbb(0:max(ks[ok]), cv[i], mu[i], rr[i])))
  }, 0)
}
snp = rbindlist(lapply(names(RULES), function(r) rbindlist(lapply(names(DS), function(tn) {
  pn = lapply(names(ARMS), function(k) arm_pnull(ARMS[[k]], rho[[k]], DS[[tn]], r))
  e = sum(1 - (1 - pn[[1]]) * (1 - pn[[2]]))
  obs = sum(pass[[paste(r, tn, sep = "__")]][ag])
  data.table(rule = r, t = tn, null = "snp_bb_ntc", ag_sites = sum(ag), ag_positives = obs, expected_chance = e,
             fdr = e / pmax(obs, 1), above_chance = obs - e)
}))))
chance = rbind(dec, snp, fill = TRUE)
setcolorder(chance, c("rule", "t", "null"))
fwrite(chance, file.path(OUT, "fisher_chance.tsv"), sep = "\t")
print(chance[, .(rule, t, null, ag_positives, expected = round(expected_chance, 1), fdr = round(fdr, 2),
                 ci = fifelse(is.na(fdr_lo), "", sprintf("%.2f-%.2f", fdr_lo, fdr_hi)), tc_rate = round(tc_rate, 3),
                 decoy_rate = round(decoy_rate, 4))])
cat("wrote", OUT, "\n")
