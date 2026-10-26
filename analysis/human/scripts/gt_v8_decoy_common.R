# =============================================================================
# gt_v8_decoy_common.R   (2026-10-17)
#
# Shared by gt_v8_empirical_null.R and gt_v8_decoy_fdr_curve.R: the v8 12-type site universe with
# per-arm KO counts, and the decoy-stratified estimate of chance A>G positives.
#   v8 universe     12-type gate sites (gt_eval_v6/spectrum12/{eval,eval_depth}/truth12_sites.tsv.gz;
#                   dbSNP156 common positions already removed) kept by the v8 cleaning: no PGP1 SNV and
#                   nearby-mismatch pass after PGP1 removal (gt_eval_v8/truth/support/truth12_site_flags.tsv.gz;
#                   gt_truth_v7_variants.R); multiallelic: A>G restricted to the v8 truth sites, other types
#                   lose positions that carry more than one mismatch type
#   groups          AG; TC (mis-stranded editing, positive control); decoy = the other 10 types
#   strata          NTC allele fraction (< 0.1, 0.1-0.2, 0.2-0.35, 0.35-0.65, 0.65-0.9, >= 0.9) x expected
#                   edited KO reads in the better-covered arm (KO cov x NTC edit: < 3, 3-5, 5-10, 10-20, >= 20);
#                   a stratum with < 30 decoys borrows the rate of its allele-fraction bin pooled over
#                   coverage, and that bin, if also < 30, the overall decoy rate
#   expected chance positives = sum over strata of n_AG(s) x decoy pass rate(s); FDR upper bound =
#                   expected / observed A>G positives (assumes every A>G site is unaffected by the KO);
#                   95% CI from bootstrap draws of the decoy sites
# Edits are rounded to 4 decimals before comparison, as classify() in gt_truth_agtc_test_v6.R.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
V8 = file.path(PM, "gt_eval_v8")
DECOY_OUT = file.path(V8, "empirical_null")
GATES = c(expected = "eval", depth = "eval_depth")
DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
MIN_DECOY = 30
LAM_BR = c(0, 3, 5, 10, 20, Inf); LAM_LAB = c("<3", "3-5", "5-10", "10-20", ">=20")
AF_BR = c(0, 0.1, 0.2, 0.35, 0.65, 0.9, Inf); AF_LAB = c("<0.1", "0.1-0.2", "0.2-0.35", "0.35-0.65", "0.65-0.9", ">=0.9")
ARMS = c(ADAR1 = "adar1", ADAR2 = "adarb1")
rnd4 = function(x) as.numeric(sprintf("%.4f", x))

load_v8_universe = function(g) {
  flags = fread(file.path(V8, "truth/support/truth12_site_flags.tsv.gz"),
                select = c("gate", "chr", "pos", "strand", "type", "pgp1", "nomm_pass_pgp1"))[gate == g]
  v8t = fread(file.path(V8, "truth/gt_v7_truth_sites.tsv.gz"), select = c("gate", "site"))[gate == g]
  d = fread(file.path(PM, "gt_eval_v6/spectrum12", GATES[[g]], "truth12_sites.tsv.gz"))
  d = merge(d, flags, by = c("chr", "pos", "strand", "type"))
  d = d[pgp1 == FALSE & nomm_pass_pgp1 == TRUE]
  d[, site := paste0(chr, ":", pos)]
  d = d[type != "AG" | site %chin% v8t$site]
  d = d[type == "AG" | !site %chin% d[, .N, by = site][N > 1, site]]
  d[, grp := fifelse(type == "AG", "AG", fifelse(type == "TC", "TC", "decoy"))]
  annotate_arms(d)
}

# per-arm KO edit, expected edited KO reads, binomial p; strata (in place)
annotate_arms = function(d) {
  d[, ntc_edit := rnd4(ntc_alt / ntc_cov)]
  for (a in ARMS) {
    cov = d[[paste0(a, "_cov")]]; alt = d[[paste0(a, "_alt")]]
    set(d, j = paste0(a, "_ke"), value = fifelse(!is.na(cov) & cov > 0, rnd4(alt / pmax(cov, 1)), NA_real_))
    set(d, j = paste0(a, "_lam"), value = fifelse(is.na(cov), 0, cov * d$ntc_edit))
    set(d, j = paste0(a, "_p"), value = fifelse(!is.na(cov) & cov > 0, pbinom(alt, cov, d$ntc_edit), 1))
  }
  d[, lam := pmax(adar1_lam, adarb1_lam)]
  d[, s_lam := cut(lam, LAM_BR, LAM_LAB, right = FALSE)]
  d[, s_af := cut(ntc_edit, AF_BR, AF_LAB, right = FALSE)]
  d[, stratum := paste(s_af, s_lam, sep = "|")]
  d[]
}

# pass: logical matrix (rows = rows of d, columns = rules / thresholds); returns one row per column
decoy_fdr = function(d, pass, n_boot = 1000, seed = 20261017) {
  pass = as.matrix(pass)
  is_ag = d$grp == "AG"; is_dc = d$grp == "decoy"
  n_ag_s = table(d$stratum[is_ag])
  obs = colSums(pass[is_ag, , drop = FALSE])
  dst = d$stratum[is_dc]; daf = as.character(d$s_af[is_dc]); P = pass[is_dc, , drop = FALSE] * 1
  sa = names(n_ag_s); af_of = sub("\\|.*", "", sa)
  est = function(idx) {
    st = dst[idx]; af = daf[idx]; X = P[idx, , drop = FALSE]
    ns = table(factor(st, levels = sa)); ks = rowsum(X, factor(st, levels = sa), reorder = TRUE)
    ks = ks[match(sa, rownames(ks)), , drop = FALSE]; ks[is.na(ks)] = 0
    na = table(factor(af, levels = AF_LAB)); ka = rowsum(X, factor(af, levels = AF_LAB), reorder = TRUE)
    ka = ka[match(AF_LAB, rownames(ka)), , drop = FALSE]; ka[is.na(ka)] = 0; rownames(ka) = AF_LAB
    overall = colMeans(X)
    r_af = sweep(ka, 1, pmax(as.numeric(na), 1), "/")
    few = as.numeric(na) < MIN_DECOY
    if (any(few)) r_af[few, ] = matrix(overall, sum(few), ncol(X), byrow = TRUE)
    r = sweep(ks, 1, pmax(as.numeric(ns), 1), "/")
    low = as.numeric(ns) < MIN_DECOY
    r[low, ] = r_af[af_of[low], , drop = FALSE]
    colSums(r * as.numeric(n_ag_s))
  }
  e = est(seq_along(dst))
  set.seed(seed)
  bt = vapply(seq_len(n_boot), function(i) est(sample.int(length(dst), replace = TRUE)), numeric(ncol(P)))
  bt = matrix(bt, nrow = ncol(P))
  data.table(rule = colnames(pass), ag_sites = sum(is_ag), decoy_sites = sum(is_dc), tc_sites = sum(d$grp == "TC"),
             ag_positives = obs, ag_rate = obs / sum(is_ag),
             tc_rate = colMeans(pass[d$grp == "TC", , drop = FALSE]), decoy_rate = colMeans(P),
             expected_chance = e, fdr = e / pmax(obs, 1),
             fdr_lo = apply(bt, 1, quantile, 0.025) / pmax(obs, 1),
             fdr_hi = apply(bt, 1, quantile, 0.975) / pmax(obs, 1),
             above_chance = obs - e)
}

strata_table = function(d) {
  d[grp %in% c("AG", "decoy"), .N, by = .(s_af, s_lam, grp)][, dcast(.SD, s_af + s_lam ~ grp, value.var = "N", fill = 0)]
}
