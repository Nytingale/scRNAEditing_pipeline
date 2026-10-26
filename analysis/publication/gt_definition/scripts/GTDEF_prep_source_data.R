#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_prep_source_data.R      (stage 1 of 2 -- ANALYSIS ONLY, writes TSVs)
#
# Ground-truth v6 thresholds.
#   step1_quality   BQ/MQ: substitution spectrum of every aligned base (merged
#                   NTC and ADAR1_KO, BQ/MQ 30 vs 20); no ground-truth threshold
#   step2_universe  dbSNP common: NTC allele fractions and ADAR1-KO response of
#                   every candidate (NTC >= 10 reads); no ground-truth threshold
#   robustness      every downstream threshold: 7-axis neighbour lattice around
#                   the centre, plus the informative-KO bound m
#   final           the frozen truth set
# 2026-09-27: the step-3 coverage diagnostics were retired (superseded by the
# robustness figure); steps 4-6 sections remain but have no sweeps and do not
# run. Earlier versions: publication/gt_definition/_superseded/scripts/.
# No ggplot here -- figures are drawn by GTDEF_gtdefinition_plot.R.
#
# INPUT directories (produced by human/pm/scripts/run_gt_truth_rebuild_v6.sh)
#   human/pm/gt_eval_v6/FROZEN_RULES.tsv, sweeps/provisional.tsv   parameter values
#   human/pm/gt_eval_v6/sweeps/<step>/grid.tsv, runs/<run_id>/      one engine run per grid row
#   human/pm/gt_eval_v6/spectrum/mismatch_spectrum.tsv             step 1: substitution counts
#   human/pm/gt_eval_v6/scan/*.hist.tsv                           counts over ALL covered A
#   human/pm/gt_eval_v6/candidates/                                cands_stats, candidates_annot
#   human/pm/gt_eval_v6/pileups/bq<BQ>/                            per-unit pileups at candidates
#   human/pm/gt_eval_v6/truth/                                     final truth set (after `final`)
#   human/pm/02_gt_bams_umi_v6/*.bam(.bai)                         idxstats only (mapped reads)
#
# OUTPUT directory
#   publication/gt_definition/report/<DATASET>/GTDEF_step<N>_*.tsv, GTDEF_final_*.tsv
#
# Usage: Rscript GTDEF_prep_source_data.R [DATASET] [STEP|all]    (default pm_gt all)
# ============================================================================

source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "GTDEF_common.R"))
# Loaded after GTDEF_common.R (which attaches data.table) so dplyr verbs win.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

argv = commandArgs(trailingOnly = TRUE)
gtdef_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
STEP_ARG = if (length(argv) >= 2) argv[2] else "all"
msg = function(...) cat("  ", ..., "\n", sep = "")
wr = function(df, stem) {
  f = file.path(REPORT_DIR, paste0(stem, ".tsv"))
  data.table::fwrite(df, f, sep = "\t")
  msg(stem, ".tsv  (", nrow(df), " rows)")
  invisible(f)
}
read_tbl = function(f, ...) as_tibble(data.table::fread(f, sep = "\t", showProgress = FALSE, ...))
want = function(st) (STEP_ARG == "all" || STEP_ARG == st) && file.exists(file.path(GT_V6_SWEEPS, st, "grid.tsv"))


# (A) Parameters and shared helpers  ----
PAR = gt_v6_params()
pv = function(k) PAR$value[[k]]
pn = function(k) as.numeric(PAR$value[[k]])
STEP_PARAMS = list(step1_quality = "bq", step2_universe = "snp_filter", step3_coverage = "min_cov",
                   step4_control_editing = c("min_edit_reads", "min_ctrl_edit"),
                   step5_reproducibility = c("min_cov_teratoma", "consistency_mode", "consistency_min"),
                   step6_depletion = "depletion")
msg("parameters: ", paste0(names(PAR$value), "=", PAR$value, ifelse(PAR$frozen, "*", ""), collapse = "  "),
    "   (* frozen)")

KO_ARMS = GT_TREAT_ARMS
SITE_COLS = c("site", "class", "powered", "depleted_ci", "in_rediportal",
              paste0("cov_", KO_ARMS), paste0("alt_", KO_ARMS))
jaccard = function(a, b) if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA_real_

read_run = function(st, run_id, cols = SITE_COLS) {
  read_tbl(file.path(GT_V6_SWEEPS, st, "runs", run_id, "gt_truth_sites.tsv.gz"), select = cols) %>%
    mutate(class = as.character(class))
}

# The step's reference run = the grid row at the current (frozen or provisional)
# values of that step's parameters; overlaps are measured against it.
ref_run_id = function(st, grid) {
  keep = rep(TRUE, nrow(grid))
  for (k in STEP_PARAMS[[st]]) keep = keep & gt_same_value(grid[[k]], pv(k))
  if (sum(keep) != 1) stop(st, ": expected exactly one reference run at current values, found ", sum(keep))
  grid$run_id[keep]
}

# One row per run: truth-set counts + overlap with the reference run.
step_metrics = function(st) {
  grid = read_tbl(file.path(GT_V6_SWEEPS, st, "grid.tsv"), colClasses = "character")
  miss = grid$run_id[!file.exists(file.path(GT_V6_SWEEPS, st, "runs", grid$run_id, "gt_truth_sites.tsv.gz"))]
  if (length(miss)) stop(st, ": runs missing: ", paste(miss, collapse = ", "))
  ref = ref_run_id(st, grid)
  runs = set_names(grid$run_id) %>% map(~ read_run(st, .x))
  ref_true = runs[[ref]] %>% filter(class == "TRUE") %>% pull(site)
  ref_pw = runs[[ref]] %>% filter(class == "TRUE", powered) %>% pull(site)
  metrics = map_dfr(grid$run_id, function(id) {
    s = runs[[id]]; tr = s %>% filter(class == "TRUE")
    mc = as.numeric(grid$min_cov[grid$run_id == id])
    # zero edited reads in every KO arm that reaches min_cov: depletion "passes"
    # at any d without evidence (the high-d plateau)
    ko_cov_ok = map(KO_ARMS, ~ tr[[paste0("cov_", .x)]] >= mc)
    ko_zero = map2(KO_ARMS, ko_cov_ok, ~ !.y | tr[[paste0("alt_", .x)]] == 0) %>% reduce(`&`)
    tibble(run_id = id,
           n_universe = nrow(s), n_true = nrow(tr),
           n_true_powered = sum(tr$powered), n_true_ci = sum(tr$depleted_ci),
           n_true_zero_ko = sum(ko_zero & reduce(ko_cov_ok, `|`)),
           n_fp_known = sum(s$class == "FP_known_not_depleted"), n_fp_novel = sum(s$class == "FP_novel"),
           n_ambiguous = sum(s$class == "ambiguous"),
           pct_true_rediportal = if (nrow(tr)) 100 * mean(tr$in_rediportal) else NA_real_,
           jaccard_vs_ref = jaccard(tr$site, ref_true),
           frac_ref_true_kept = if (length(ref_true)) mean(ref_true %in% tr$site) else NA_real_,
           jaccard_powered_vs_ref = jaccard(tr$site[tr$powered], ref_pw))
  })
  grid %>% left_join(metrics, by = "run_id") %>%
    mutate(n_fp = n_fp_known + n_fp_novel,
           is_ref = run_id == ref,
           status = if (all(PAR$frozen[STEP_PARAMS[[st]]])) "adopted" else "provisional")
}

# Genome-wide control yield from the scan histograms (every covered A).
hist_all = function() {
  list.files(GT_V6_SCAN, "\\.hist\\.tsv$", full.names = TRUE) %>% map_dfr(read_tbl) %>%
    group_by(cov_bin_min, edited_reads_capped10, edit_bin_min) %>%
    summarise(n_positions = sum(n_positions), .groups = "drop")
}
hist_pass = function(h, cov = 1, k = 0, e = 0)
  sum(h$n_positions[h$cov_bin_min >= cov & h$edited_reads_capped10 >= k & h$edit_bin_min >= e - 1e-9])

read_pile = function(unit, bq = pv("bq")) {
  read_tbl(file.path(GT_V6_PILEUPS, paste0("bq", bq), paste0("pileup_", unit, ".tsv")),
           select = c("coverage", "variant_cov", "editing_level")) %>%
    transmute(cov = coverage, alt = variant_cov, edit = as.numeric(editing_level))
}


# (Step 1) Pileup quality: substitution spectrum of every aligned base  ----
# No ground-truth threshold is used. Input: spectrum/mismatch_spectrum.tsv
# (run_gt_truth_rebuild_v6.sh spectrum*): strand-aware base counts over every
# aligned base of merged NTC and merged ADAR1_KO at BQ/MQ 30 and 20.
SPEC_F = file.path(GT_V6_DIR, "spectrum", "mismatch_spectrum.tsv")
if (STEP_ARG %in% c("all", "step1_quality") && file.exists(SPEC_F)) {
  msg("step1_quality (substitution spectrum)")
  sp = read_tbl(SPEC_F, integer64 = "double")   # genome-wide counts exceed 2^31
  # rate of each RNA-strand substitution per 10^6 aligned bases of that reference base
  rates = sp %>% group_by(sample, bq_mq, rna_ref) %>% mutate(n_ref_bases = sum(n)) %>% ungroup() %>%
    filter(rna_ref != rna_base) %>%
    mutate(type = paste0(rna_ref, ">", rna_base), per_million = 1e6 * n / n_ref_bases) %>%
    group_by(sample, bq_mq) %>% mutate(pct_of_mismatches = 100 * n / sum(n)) %>% ungroup()
  wr(rates %>% select(sample, bq_mq, type, rna_ref, rna_base, n, n_ref_bases, per_million, pct_of_mismatches), "GTDEF_step1_spectrum")
  # ADAR dependence of each substitution: ADAR1_KO rate / NTC rate at the same quality
  wr(rates %>% select(sample, bq_mq, type, per_million) %>%
       pivot_wider(names_from = sample, values_from = per_million) %>%
       mutate(ko_over_ntc = merged_ADAR1_KO / merged_NTC_ctrl),
     "GTDEF_step1_ntc_vs_ko")
  wr(sp %>% group_by(sample, bq_mq) %>%
       summarise(aligned_bases = sum(n), mismatches = sum(n[rna_ref != rna_base]),
                 a_to_g = sum(n[rna_ref == "A" & rna_base == "G"]), .groups = "drop") %>%
       mutate(a_to_g_share = a_to_g / mismatches,
              non_ag_per_million = 1e6 * (mismatches - a_to_g) / aligned_bases),
     "GTDEF_step1_summary")
}


# (Step 2) Site universe: dbSNP common sites vs the rest  ----
# No ground-truth threshold is used: every candidate (the upstream universe:
# NTC >= 5 reads, >= 1 edited read, >= 1% editing) with a display floor of
# NTC >= 10 reads so allele fractions are stable.
STEP2_NTC_MIN = 10
if (STEP_ARG %in% c("all", "step2_universe")) {
  msg("step2_universe (allele fractions and ADAR1-KO response)")
  ann = read_tbl(file.path(GT_V6_DIR, "candidates", "candidates_annot.tsv")) %>%
    mutate(group = if_else(in_dbsnp_common == 1, "Sites in dbSNP", "Sites not in dbSNP"))   # dbSNP156 common
  pile = function(unit) read_tbl(file.path(GT_V6_PILEUPS, "bq30", paste0("pileup_", unit, ".tsv")),
                                 select = c("chr", "pos", "coverage", "variant_cov")) %>%
    transmute(site = paste0(chr, ":", pos), cov = coverage, alt = variant_cov)
  s = ann %>% select(site, group) %>%
    inner_join(pile("merged_NTC_ctrl") %>% rename(ntc_cov = cov, ntc_alt = alt), by = "site") %>%
    left_join(pile("merged_ADAR1_KO") %>% rename(ko_cov = cov, ko_alt = alt), by = "site") %>%
    mutate(across(c(ko_cov, ko_alt), ~ coalesce(.x, 0)), ntc_af = ntc_alt / ntc_cov) %>%
    filter(ntc_cov >= STEP2_NTC_MIN)
  wr(s %>% group_by(group) %>%
       summarise(n_sites = n(),
                 median_af = median(ntc_af), pct_af_ge_0.35 = 100 * mean(ntc_af >= 0.35), .groups = "drop") %>%
       mutate(ntc_min_reads = STEP2_NTC_MIN),
     "GTDEF_step2_summary")
  # (a) NTC allele-fraction distribution, share of sites per 2% bin within each group
  wr(s %>% mutate(bin = pmin(floor(ntc_af / 0.02), 49) * 0.02 + 0.01) %>%
       count(group, bin, name = "n") %>% group_by(group) %>% mutate(frac = n / sum(n)) %>% ungroup(),
     "GTDEF_step2_af_hist")
  # (b) ADAR1-KO response: pooled KO editing vs pooled NTC editing in NTC-editing
  # bins (pooling reads avoids any per-site KO coverage floor); Wilson 95% CI
  BR = c(0.01, 0.05, 0.10, 0.20, 0.35, 0.65, 0.90, 1.0001)
  wr(s %>% mutate(af_bin = cut(ntc_af, BR, right = FALSE)) %>% filter(!is.na(af_bin)) %>%
       group_by(group, af_bin) %>%
       summarise(n_sites = n(), ntc_edit = sum(ntc_alt) / sum(ntc_cov),
                 ko_alt = sum(ko_alt), ko_cov = sum(ko_cov), .groups = "drop") %>%
       mutate(ko_edit = if_else(ko_cov > 0, ko_alt / ko_cov, NA_real_),
              z = qnorm(0.975), den = 1 + z^2 / ko_cov,
              ko_lo = (ko_edit + z^2 / (2 * ko_cov) - z * sqrt(ko_edit * (1 - ko_edit) / ko_cov + z^2 / (4 * ko_cov^2))) / den,
              ko_hi = (ko_edit + z^2 / (2 * ko_cov) + z * sqrt(ko_edit * (1 - ko_edit) / ko_cov + z^2 / (4 * ko_cov^2))) / den,
              ko_over_ntc = ko_edit / ntc_edit) %>%
       select(-z, -den),
     "GTDEF_step2_ko_response")
}


# (Step 4) Control editing: min_edit_reads x min_ctrl_edit  ----
if (want("step4_control_editing")) {
  msg("step4_control_editing")
  m = step_metrics("step4_control_editing") %>%
    mutate(k = as.integer(min_edit_reads), pct = as.numeric(min_ctrl_edit))
  wr(m, "GTDEF_step4_control_editing")
  # genome-wide: covered A passing the control criteria at the chosen min_cov
  h = hist_all()
  wr(expand_grid(k = 1:5, pct = c(0.01, 0.02, 0.05, 0.10)) %>%
       mutate(min_cov = pn("min_cov"),
              n_positions = map2_dbl(k, pct, ~ hist_pass(h, cov = pn("min_cov"), k = .x, e = .y))),
     "GTDEF_step4_ctrl_yield")
  # control editing distribution at candidates, by coverage stratum
  ctl = read_pile(paste0("merged_", GT_CTRL_ARM)) %>% filter(alt > 0)
  strata = sort(unique(c(pn("min_cov"), 20, 50, 100)))
  wr(map_dfr(strata, function(cv) {
       d = ctl %>% filter(cov >= cv) %>% arrange(edit) %>% mutate(cum = row_number() / n())
       d %>% slice(unique(round(seq(1, n(), length.out = min(n(), 3000))))) %>% transmute(cov_min = cv, edit, cum)
     }), "GTDEF_step4_editing_ecdf")
  wr(expand_grid(cov_min = strata, k = 1:5) %>% mutate(floor_edit = k / cov_min), "GTDEF_step4_editing_floors")
  n_tests = m$n_universe[m$is_ref]
  wr(tibble(coverage = unique(round(exp(seq(log(pn("min_cov")), log(1000), length.out = 120))))) %>%
       mutate(k_binom_uncorrected = qbinom(1 - GT_ALPHA, coverage, GT_SEQ_ERROR) + 1L,
              k_binom_bonferroni = qbinom(1 - GT_ALPHA / n_tests, coverage, GT_SEQ_ERROR) + 1L,
              n_tests = n_tests),
     "GTDEF_step4_binomial_requirement")
}


# (Step 5) Reproducibility: min_cov_teratoma x consistency rule  ----
if (want("step5_reproducibility")) {
  msg("step5_reproducibility")
  m = step_metrics("step5_reproducibility") %>%
    mutate(tcov = as.integer(min_cov_teratoma),
           rule = if_else(consistency_mode == "evaluable", "evaluable", paste0("strict_ge", consistency_min)))
  wr(m, "GTDEF_step5_reproducibility")
  ref = m$run_id[m$is_ref]
  # UMI depth per unit (index only)
  idx_reads = function(bam) {
    x = system2("samtools", c("idxstats", shQuote(bam)), stdout = TRUE)
    sum(as.numeric(vapply(strsplit(x, "\t"), `[`, "", 3)))
  }
  depth = tibble(unit = sub("\\.bam$", "", list.files(GT_V6_BAMDIR, "\\.bam$"))) %>%
    mutate(mapped_reads = map_dbl(unit, ~ idx_reads(file.path(GT_V6_BAMDIR, paste0(.x, ".bam")))),
           tier = if_else(startsWith(unit, "merged_"), "merged", "sample"),
           teratoma = if_else(tier == "sample", sub("_.*$", "", unit), NA_character_),
           arm = sub("^(merged|T[0-9]+)_", "", unit))
  wr(depth, "GTDEF_step5_unit_depth")
  # evaluability / positivity per teratoma, at the reference teratoma depth
  full = read_tbl(file.path(GT_V6_SWEEPS, "step5_reproducibility", "runs", ref, "gt_truth_sites.tsv.gz"))
  ev = full %>% filter(.data[[paste0("cov_", GT_CTRL_ARM)]] >= pn("min_cov"))
  wr(map_dfr(GT_TERATOMAS, ~ tibble(teratoma = .x, n_sites = nrow(ev),
                                    n_evaluable = sum(ev[[paste0("eval_", .x)]] == 1L),
                                    n_positive = sum(ev[[paste0("pos_", .x)]] == 1L))) %>%
       left_join(depth %>% filter(tier == "sample", arm == GT_CTRL_ARM) %>% select(teratoma, mapped_reads), by = "teratoma") %>%
       mutate(evaluable_rate = n_evaluable / n_sites,
              positive_rate_of_evaluable = if_else(n_evaluable > 0, n_positive / n_evaluable, NA_real_),
              min_cov_teratoma = pv("min_cov_teratoma")),
     "GTDEF_step5_teratoma_positivity")
  wr(ev %>% count(n_evaluable, name = "n_sites") %>% mutate(min_cov_teratoma = pv("min_cov_teratoma")),
     "GTDEF_step5_reachability")
  # rule agreement on powered TRUE sets, at the reference teratoma depth
  at_t = m %>% filter(gt_same_value(min_cov_teratoma, pv("min_cov_teratoma")))
  sets = set_names(at_t$run_id, at_t$rule) %>%
    map(~ read_run("step5_reproducibility", .x) %>% filter(class == "TRUE", powered) %>% pull(site))
  wr(expand_grid(rule_a = names(sets), rule_b = names(sets)) %>%
       mutate(n_a = map_int(rule_a, ~ length(sets[[.x]])), n_b = map_int(rule_b, ~ length(sets[[.x]])),
              n_shared = map2_int(rule_a, rule_b, ~ length(intersect(sets[[.x]], sets[[.y]]))),
              jaccard = map2_dbl(rule_a, rule_b, ~ jaccard(sets[[.x]], sets[[.y]])),
              frac_of_smaller = if_else(pmin(n_a, n_b) > 0, n_shared / pmin(n_a, n_b), NA_real_),
              min_cov_teratoma = pv("min_cov_teratoma")),
     "GTDEF_step5_rule_agreement")
}


# (Step 6) Depletion d, with statistical power  ----
if (want("step6_depletion")) {
  msg("step6_depletion")
  m = step_metrics("step6_depletion") %>% mutate(d = as.numeric(depletion))
  wr(m, "GTDEF_step6_depletion")
  ref = m$run_id[m$is_ref]
  full = read_tbl(file.path(GT_V6_SWEEPS, "step6_depletion", "runs", ref, "gt_truth_sites.tsv.gz"))
  gated = full %>% filter(evaluable, ctrl_edited, consistent)   # independent of d
  wr(map_dfr(c(GT_CTRL_ARM, KO_ARMS), function(arm) {
       g = gated[[paste0("cov_", arm)]]
       tibble(arm = arm, median_cov_gated = median(g), q25_gated = quantile(g, .25),
              q75_gated = quantile(g, .75), n_gated = length(g))
     }), "GTDEF_step6_reads_available")
  reads_needed = function(target) if (!is.finite(target) || target <= 0) NA_integer_ else
    as.integer(ceiling(log(GT_ALPHA) / log(1 - target)))
  wr(expand_grid(ctrl_edit = seq(0.1, 0.9, 0.1), depletion = GT_DEP_SERIES) %>%
       mutate(target = (1 - depletion) * ctrl_edit, reads_needed = map_int(target, reads_needed)),
     "GTDEF_step6_power_reads")
  # largest decrease each gated site's own KO depth could demonstrate
  # (zero edited KO reads: 95% one-sided CP upper bound = 1 - alpha^(1/n));
  # "powered at d" <=> d_max >= d
  cp0 = function(n) ifelse(n > 0, 1 - GT_ALPHA^(1 / n), 1)
  d_arm = function(n, ctrl) ifelse(n >= pn("min_cov") & ctrl > 0, pmax(0, 1 - cp0(n) / ctrl), NA_real_)
  dm = gated %>% transmute(ctrl = .data[[paste0("edit_", GT_CTRL_ARM)]],
                           a1 = d_arm(.data[[paste0("cov_", KO_ARMS[1])]], ctrl),
                           a2 = d_arm(.data[[paste0("cov_", KO_ARMS[2])]], ctrl)) %>%
    mutate(either = pmax(coalesce(a1, -1), coalesce(a2, -1)))
  wr(tibble(d = seq(0, 1, 0.01)) %>%
       mutate(n_gated = nrow(dm),
              frac_either = map_dbl(d, ~ mean(dm$either >= .x)),
              frac_arm1 = map_dbl(d, ~ mean(!is.na(dm$a1) & dm$a1 >= .x)),
              frac_arm2 = map_dbl(d, ~ mean(!is.na(dm$a2) & dm$a2 >= .x)),
              arm1 = KO_ARMS[1], arm2 = KO_ARMS[2]),
     "GTDEF_step6_demonstrable")
}


# (Robustness) 7-axis neighbour lattice around the fixed centre  ----
# Inputs from human/scripts/gt_robustness_grid_v6.R in human/pm/gt_eval_v6/robustness/
# Observed Jaccard values only; no stability cut-off is imposed (user, 2026-09-27).
RB = file.path(GT_V6_DIR, "robustness")
if (STEP_ARG %in% c("all", "robustness") && file.exists(file.path(RB, "combos.tsv.gz"))) {
  msg("robustness")
  # Per-teratoma NTC coverage is an axis (8-12), swept like NTC coverage.
  cb = read_tbl(file.path(RB, "combos.tsv.gz"))
  nball = read_tbl(file.path(RB, "neighbours.tsv.gz")) %>% filter(combo %in% cb$combo)
  # continuous thresholds vs the discrete ones (edited reads: integer steps of
  # +-33%; rule: a whole replicate)
  CONT = c("min_cov", "min_cov_teratoma", "ko_min_expected", "min_ctrl_edit", "depletion")
  agg = nball %>% group_by(combo) %>%
    summarise(min_jaccard_1 = min(jaccard[abs(step) == 1]), min_jaccard_2 = min(jaccard),
              min_jaccard_1_cont = min(jaccard[abs(step) == 1 & axis %in% CONT]),
              min_jaccard_1_disc = min(jaccard[abs(step) == 1 & !axis %in% CONT]),
              n_neighbours_1 = sum(abs(step) == 1), .groups = "drop")
  cb = cb %>% select(-any_of(c("min_jaccard_1", "min_jaccard_2", "n_neighbours_1", "stable_1", "stable_2",
                                "interior_1", "interior_2", "max_frac_lost_1", "max_frac_lost_2",
                                "min_neighbour_sites_1", "n_neighbours_2"))) %>%
    left_join(agg, by = "combo")
  wr(cb, "GTDEF_robust_combos")
  centre = cb %>% filter(is_centre)
  nb = nball %>% filter(combo == centre$combo)
  axis_centre = c(min_cov = "10", min_cov_teratoma = "10", ko_min_expected = "3", min_edit_reads = "3", min_ctrl_edit = "0.05",
                  depletion = "0.5", rule = "evaluable")
  wr(bind_rows(nb,
               tibble(combo = centre$combo, axis = names(axis_centre), step = 0L, neighbour = centre$combo,
                      neighbour_value = unname(axis_centre), n_true = centre$n_true, n_neighbour = centre$n_true,
                      jaccard = 1, frac_lost = 0)) %>%
       arrange(axis, step),
     "GTDEF_robust_centre_neighbours")
  if (file.exists(file.path(RB, "parity.tsv"))) wr(read_tbl(file.path(RB, "parity.tsv")), "GTDEF_robust_parity")
}


# (Robustness, m) What the informative-KO threshold m means  ----
# m = expected edited KO reads if the KO did NOT deplete editing = KO coverage n
# x merged-NTC editing p. With no depletion each KO read is edited with
# probability p, so P(0 edited KO reads) = (1 - p)^n <= exp(-n p) = exp(-m);
# m = -log(alpha) bounds that chance by alpha (0.05 -> 3.00, 0.01 -> 4.61).
# Theory tables need no data; the observed curve is robustness/m_curve.tsv
# (gt_robustness_grid_v6.R mcurve, every other threshold at the centre).
if (STEP_ARG %in% c("all", "robustness")) {
  msg("robustness: informative-KO m")
  P_SHOW = c(0.05, 0.10, 0.20, 0.50)
  wr(tidyr::expand_grid(m = seq(0.05, 6, 0.05), ntc_edit = P_SHOW) %>%
       mutate(ko_cov_min = ceiling(m / ntc_edit - 1e-9),
              p0_exact = (1 - ntc_edit)^ko_cov_min,
              p0_bound = exp(-m)),
     "GTDEF_robust_m_theory")
  wr(tidyr::expand_grid(ntc_edit = exp(seq(log(0.03), log(0.8), length.out = 120)), alpha = c(0.05, 0.01)) %>%
       mutate(m = -log(alpha),
              ko_cov_m_rule = ceiling(m / ntc_edit - 1e-9),
              ko_cov_exact = ceiling(log(alpha) / log(1 - ntc_edit) - 1e-9)),
     "GTDEF_robust_ko_cov_needed")
  mf = file.path(RB, "m_curve.tsv")
  if (file.exists(mf)) wr(read_tbl(mf), "GTDEF_robust_m_curve")
}


# (v8) Robustness and benchmark sweeps on the v8 site set  ----
# Robustness: human/scripts/gt_robustness_v8.R (A>G read orientation; multiallelic, PGP1 and
# nearby-mismatch positions removed; both gates). Benchmark: gt_eval_v8/compare_config_rank.tsv,
# call state mismatch, reported configurations (F1 per tier = mean over units, gt_score_v7.R).
RB8 = file.path(V8_DIR, "robustness")
if (STEP_ARG %in% c("all", "v8") && file.exists(file.path(RB8, "centre_neighbours.tsv"))) {
  msg("v8")
  par8 = read_tbl(file.path(RB8, "parity.tsv"))
  if (!all(par8$identical)) stop("gt_robustness_v8.R centre does not reproduce the v8 positives")
  wr(read_tbl(file.path(RB8, "centre_neighbours.tsv")), "GTDEF_robust_centre_neighbours_v8")
  wr(read_tbl(file.path(RB8, "m_curve.tsv")), "GTDEF_robust_m_curve_v8")
  wr(read_tbl(file.path(RB8, "ko_cov_curve.tsv")), "GTDEF_robust_ko_cov_curve_v8")
  wr(read_tbl(file.path(RB8, "cleaning_counts.tsv")), "GTDEF_robust_cleaning_v8")
  wr(par8, "GTDEF_robust_parity_v8")
  cr = read_tbl(V8_CMP_RANK) %>% filter(state == V8_STATE) %>%
    select(gate, t, config, caller, f1_arm, f1_sample, f1_celltype) %>%
    pivot_longer(starts_with("f1_"), names_to = "tier", names_prefix = "f1_", values_to = "f1")
  wr(cr %>% group_by(gate, tier, t, caller) %>% slice_max(f1, n = 1, with_ties = FALSE) %>% ungroup() %>%
       select(gate, tier, t, caller, config, f1),
     "GTDEF_sweep_caller_best_v8")
  ts = sort(unique(cr$t))
  wr(map_dfr(combn(ts, 2, simplify = FALSE), function(p)
       cr %>% filter(t %in% p) %>% pivot_wider(names_from = t, values_from = f1) %>%
         group_by(gate, tier) %>%
         summarise(t_a = p[1], t_b = p[2], n_config = n(),
                   rho_f1 = cor(.data[[p[1]]], .data[[p[2]]], method = "spearman", use = "complete.obs"),
                   .groups = "drop")),
     "GTDEF_sweep_rank_stability_v8")
}


# (Final) The frozen truth set  ----
tf = file.path(GT_V6_DIR, "truth", "gt_truth_sites.tsv.gz")
if ((STEP_ARG %in% c("all", "final")) && file.exists(tf)) {
  msg("final")
  h = hist_all(); n_cov = hist_pass(h)
  cs = read_tbl(file.path(GT_V6_DIR, "candidates", "cands_stats.tsv"))
  cls = read_tbl(tf, select = "class") %>% mutate(class = as.character(class)) %>% count(class, name = "N")
  wr(bind_rows(cls,
               if (pv("snp_filter") == "excluded")
                 tibble(class = "dbSNP_common_excluded", N = cs$n[cs$metric == "in_dbsnp_common"]),
               tibble(class = "below_candidate_floor", N = n_cov - cs$n[cs$metric == "candidates_all"])) %>%
       mutate(n_covered_A = n_cov),
     "GTDEF_final_truth_classes")
}

cat("\nprep complete (", STEP_ARG, ") -> ", REPORT_DIR, "\n", sep = "")
