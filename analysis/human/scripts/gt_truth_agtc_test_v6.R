#!/usr/bin/env Rscript
# =============================================================================
# gt_truth_agtc_test_v6.R   (2026-10-06)
#
# Test of strand-agnostic ground-truth definitions (HANDOFF 2026-10-05 evening, section 4-5;
# GTDEF_panel_agtc_truth): T>C positives are gene-strand A-to-I seen in reversed-orientation
# reads. The frozen v6 rules (FROZEN_RULES.tsv, gt_spectrum12_eval_v6.R) are applied to
#   AG      RNA-strand A>G only (current v6 truth; must reproduce truth_sweep_classes exactly)
#   AGTC    A>G or T>C per read orientation, keyed on chr:pos; where both orientations are
#           candidates, the one with more merged-NTC reads is kept (ties: A>G)
#   pooled  genomic A:T>G:C transition with reads of both orientations summed per unit:
#           ref A: cov = '.' + ',' + G + g, alt = G + g;  ref T: cov = '.' + ',' + C + c, alt = C + c
# each with and without removal of multiallelic sites (merged NTC: a third allele, not ref and not
# the edited base, with >= 3 reads and >= 5% of depth = cov + third-allele reads; on the site's read
# orientation for AG / AGTC, both orientations for pooled; gt_multiallelic_counts_v6.R rule).
# Removed sites are unscored (not eligible). Gates: expected-gated and depth-gated (s12).
#
# Scoring (published call sets, both window filters, *_pre100.tsv.gz; reported configurations):
#   call match AG    calls with varType AG (current)
#   call match AGTC  calls with varType AG or TC (strand assignment ignored)
#   TP / FP on positives(t) / informative negatives by chr:pos; precision = TP / (TP + FP)
#   (NA if 0), recall = TP / n positives; tier value = mean over units; F1 = 2PR / (P + R) of
#   the tier values; caller rank by cross-tier rank sum (TERATOMA_window_state_ranksum_plot.R).
# Parity: truth AG x call match AG x multiallelic kept = published gt_tier_scoring.tsv.
#
# Inputs : gt_eval_v6/spectrum12/{positions_noSNP.bed, candidates_all.tsv.gz, pileups/}
#          gt_eval_v6/scoring/sweep_infneg{,_depth}/truth_sweep_classes.tsv.gz and
#            teratoma_inputs/invariance/d0.{25,50,75}/gt_tier_scoring.tsv
#          human/pm/results_gt{,_persample,_celltype}_v6/<NTC unit>/04_filtered/**/*_pre100.tsv.gz
# Outputs: gt_eval_v6/agtc_test/
#          truth_sites.tsv.gz     variant x gate x site: type, orientation, classes, multiallelic
#          truth_counts.tsv       class sizes per variant x gate x multiallelic handling
#          tier_scores.tsv        scenario x truth x tier x configuration
#          caller_rank.tsv        caller rank per scenario x truth
#          strand_mode_recall.tsv median recall per caller x strand mode
#          parity.tsv
# Usage  : qsub human/pm/scripts/gt_truth_agtc_test_v6.sh
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
V6 = file.path(PM, "gt_eval_v6")
SP = Sys.getenv("GT_SP12_DIR", file.path(V6, "spectrum12"))
OUT = Sys.getenv("GT_AGTC_DIR", file.path(V6, "agtc_test")); dir.create(OUT, showWarnings = FALSE)
# GT_TRUTH_ONLY=1 (split-by-cells test): write truth_sites + truth_counts only; no parity with the full-data
# sweep classes (a different cell set) and no call scoring
TRUTH_ONLY = Sys.getenv("GT_TRUTH_ONLY") == "1"
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

fr = fread(file.path(V6, "FROZEN_RULES.tsv"))
P = setNames(fr$value, fr$parameter)
MIN_COV = as.numeric(P[["min_cov"]]); MIN_COV_KO = as.numeric(P[["min_cov_ko"]])
M_EXP = as.numeric(P[["ko_min_expected"]]); K = as.numeric(P[["min_edit_reads"]])
E = as.numeric(P[["min_ctrl_edit"]]); TCOV = as.numeric(P[["min_cov_teratoma"]])
stopifnot(P[["consistency_mode"]] == "evaluable", P[["snp_filter"]] == "excluded", P[["bq"]] == "30")
MIN_COV_KO_DEPTH = 10
NEG_MAX_DEP = 0.25; CI_ALPHA = 0.05
MIN_THIRD_READS = 3; MIN_THIRD_FRAC = 0.05
DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
GATES = c(expected = "sweep_infneg", depth = "sweep_infneg_depth")
# GT_GATES=<comma list> (2026-10-17): gates to build; depth<N> = the depth rule with KO cov >= N (needs GT_TRUTH_ONLY=1)
if (nzchar(Sys.getenv("GT_GATES"))) {
  gl = strsplit(Sys.getenv("GT_GATES"), ",", fixed = TRUE)[[1]]
  stopifnot(all(gl == "expected" | grepl("^depth[0-9]*$", gl)))
  GATES = setNames(ifelse(gl == "expected", "sweep_infneg", ifelse(gl == "depth", "sweep_infneg_depth", NA)), gl)
  stopifnot(TRUTH_ONLY || !anyNA(GATES))
}
ko_floor = function(gate_type) if (gate_type == "depth") MIN_COV_KO_DEPTH else as.numeric(sub("^depth", "", gate_type))
KO = c("ADAR1_KO", "ADARB1_KO"); TT = c("T1", "T2", "T3", "T4")
UNITS = c("merged_NTC_ctrl", paste0("merged_", KO), paste0(TT, "_NTC_ctrl"))
COMP = c(A = "T", C = "G", G = "C", T = "A")
BASES = c("A", "C", "G", "T")

# ---- pileups (as gt_spectrum12_eval_v6.R) ------------------------------------------------------
pos = fread(file.path(SP, "positions_noSNP.bed"), header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
pos[, idx := .I]
setkey(pos, chr, pos)
N = nrow(pos); msg("positions: ", N)
cols = c("ref_plus", "ref_minus", "A", "C", "G", "T", "a", "c", "g", "t")
read_unit = function(u) {
  fs = list.files(file.path(SP, "pileups"), pattern = paste0("^", u, "__.*\\.tsv\\.gz$"), full.names = TRUE)
  d = rbindlist(lapply(fs, function(f) fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", showProgress = FALSE)))
  d = pos[d, on = .(chr, pos), nomatch = NULL]
  m = matrix(0L, N, length(cols), dimnames = list(NULL, cols))
  for (cn in cols) m[d$idx, cn] = as.integer(d[[cn]])
  ref = rep(NA_character_, N); ref[d$idx] = d$ref
  msg(u, ": ", nrow(d), " covered positions from ", length(fs), " chunks")
  list(m = m, ref = ref)
}
U = setNames(lapply(UNITS, read_unit), UNITS)
REF = U[["merged_NTC_ctrl"]]$ref
for (u in UNITS) REF[is.na(REF)] = U[[u]]$ref[is.na(REF)]
REF = toupper(REF)
rnd4 = function(x) as.numeric(sprintf("%.4f", x))
mk = function(cov, alt) list(cov = cov, alt = alt, edit = ifelse(cov > 0, rnd4(alt / pmax(cov, 1)), 0))

# per read orientation: genomic ref / alt columns on that orientation's reads
orient_counts = function(u, idx, strand, R, Y) {
  m = U[[u]]$m; plus = strand == "+"
  alt = ifelse(plus, m[cbind(idx, match(Y, cols))], m[cbind(idx, match(tolower(COMP[[Y]]), cols))])
  refn = ifelse(plus, m[idx, "ref_plus"], m[idx, "ref_minus"])
  mk(refn + alt, alt)
}
orient_third = function(idx, strand, R, Y) {
  m = U[["merged_NTC_ctrl"]]$m; plus = strand == "+"
  th = setdiff(BASES, c(R, Y))
  t1 = ifelse(plus, m[cbind(idx, match(th[1], cols))], m[cbind(idx, match(tolower(COMP[[th[1]]]), cols))])
  t2 = ifelse(plus, m[cbind(idx, match(th[2], cols))], m[cbind(idx, match(tolower(COMP[[th[2]]]), cols))])
  list(t1 = t1, t2 = t2)
}
# pooled over both orientations at genomic ref A (alt G) or T (alt C)
pooled_counts = function(u, idx, gref) {
  m = U[[u]]$m; isA = gref == "A"
  alt = ifelse(isA, m[idx, "G"] + m[idx, "g"], m[idx, "C"] + m[idx, "c"])
  mk(m[idx, "ref_plus"] + m[idx, "ref_minus"] + alt, alt)
}
pooled_third = function(idx, gref) {
  m = U[["merged_NTC_ctrl"]]$m; isA = gref == "A"
  list(t1 = ifelse(isA, m[idx, "C"] + m[idx, "c"], m[idx, "A"] + m[idx, "a"]),
       t2 = ifelse(isA, m[idx, "T"] + m[idx, "t"], m[idx, "G"] + m[idx, "g"]))
}

cp_lower = function(k, n, a = CI_ALPHA) { o = rep(0, length(k)); ok = n > 0 & k > 0; o[ok] = qbeta(a, k[ok], n[ok] - k[ok] + 1); o }
classify = function(S, gate_type) {
  ctrl = S[["merged_NTC_ctrl"]]
  cov_arm = lapply(KO, function(a) { k = S[[paste0("merged_", a)]]
    if (gate_type == "expected") ctrl$cov >= MIN_COV & k$cov >= MIN_COV_KO & k$cov * ctrl$edit >= M_EXP
    else ctrl$cov >= MIN_COV & k$cov >= ko_floor(gate_type) })
  ctrl_ed = ctrl$alt >= K & ctrl$edit >= E
  n_eval = Reduce(`+`, lapply(TT, function(t) S[[paste0(t, "_NTC_ctrl")]]$cov >= TCOV))
  n_pos = Reduce(`+`, lapply(TT, function(t) { x = S[[paste0(t, "_NTC_ctrl")]]; x$cov >= TCOV & x$alt >= K & x$edit >= E }))
  gate = Reduce(`|`, cov_arm) & ctrl_ed & n_pos >= 2 & n_pos >= pmax(2, n_eval - 1)
  nd = lapply(seq_along(KO), function(i) { k = S[[paste0("merged_", KO[i])]]
    cov_arm[[i]] & cp_lower(k$alt, k$cov) > (1 - NEG_MAX_DEP) * ctrl$edit })
  o = data.table(eligible = gate, informative_negative = gate & nd[[1]] & nd[[2]], ntc_cov = ctrl$cov)
  for (nm in names(DS)) {
    target = (1 - DS[[nm]]) * ctrl$edit
    dep = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & S[[paste0("merged_", KO[i])]]$edit <= target))
    set(o, j = paste0("pos_", nm), value = gate & dep)
  }
  stopifnot(!any(o$informative_negative & (o$pos_t25 | o$pos_t50 | o$pos_t75)))
  o
}
is_multi = function(th, cov) { mx = pmax(th$t1, th$t2); mx >= MIN_THIRD_READS & mx >= MIN_THIRD_FRAC * (cov + th$t1 + th$t2) }

# ---- truth variants ------------------------------------------------------------------------------
cand = fread(cmd = paste("gzip -cd", shQuote(file.path(SP, "candidates_all.tsv.gz"))), sep = "\t", header = FALSE,
             select = 1:5, col.names = c("chr", "pos", "strand", "rna_ref", "rna_alt"), showProgress = FALSE)
cand = pos[cand, on = .(chr, pos), nomatch = NULL]
cand[, type := paste0(rna_ref, rna_alt)]
cand = unique(cand[type %in% c("AG", "TC"), .(chr, pos, idx, strand, type)])
msg("AG / TC candidate (position, orientation) rows: ", paste(cand[, .N, by = type][, paste(type, N)], collapse = ", "))

ori = rbindlist(lapply(c("AG", "TC"), function(ty) {
  R = substr(ty, 1, 1); Y = substr(ty, 2, 2); cd = cand[type == ty]
  gref = ifelse(cd$strand == "+", R, COMP[[R]])
  stopifnot(all(REF[cd$idx] == gref, na.rm = TRUE))
  S = lapply(setNames(UNITS, UNITS), function(u) orient_counts(u, cd$idx, cd$strand, R, Y))
  multi = is_multi(orient_third(cd$idx, cd$strand, R, Y), S[["merged_NTC_ctrl"]]$cov)
  rbindlist(lapply(names(GATES), function(g) cbind(cd[, .(site = paste0(chr, ":", pos), type, strand)],
                                                  gate = g, multiallelic = multi, classify(S, g))))
}))
ori = ori[eligible == TRUE]

pidx = which(REF %in% c("A", "T"))
Sp = lapply(setNames(UNITS, UNITS), function(u) pooled_counts(u, pidx, REF[pidx]))
pmulti = is_multi(pooled_third(pidx, REF[pidx]), Sp[["merged_NTC_ctrl"]]$cov)
pool = rbindlist(lapply(names(GATES), function(g)
  data.table(site = paste0(pos$chr[match(pidx, pos$idx)], ":", pos$pos[match(pidx, pos$idx)]),
             type = fifelse(REF[pidx] == "A", "A:G", "T:C"), strand = ".", gate = g, multiallelic = pmulti,
             classify(Sp, g))))
pool = pool[eligible == TRUE]
rm(Sp); invisible(gc())

agtc = copy(ori)[order(gate, site, -ntc_cov, type)][, .SD[1], by = .(gate, site)]
truth = rbindlist(list(ori[type == "AG"][, variant := "AG"], agtc[, variant := "AGTC"], pool[, variant := "pooled"]),
                  use.names = TRUE)
truth[, ntc_cov := NULL]
fwrite(truth, file.path(OUT, "truth_sites.tsv.gz"), sep = "\t")

# parity: AG must equal the frozen sweep classes
if (!TRUTH_ONLY) {
par = rbindlist(lapply(names(GATES), function(g) {
  sw = fread(file.path(V6, "scoring", GATES[[g]], "truth_sweep_classes.tsv.gz"))[eligible == TRUE]
  me = truth[variant == "AG" & gate == g]
  rbindlist(lapply(c("eligible", "informative_negative", paste0("pos_", names(DS))), function(cl) {
    ref = if (cl == "eligible") sw$site else if (cl == "informative_negative") sw[informative_negative == TRUE, site]
          else sw[get(sub("^pos_", "class_", cl)) == "positive", site]
    mine = if (cl == "eligible") me$site else me[get(cl) == TRUE, site]
    data.table(check = "truth_classes", gate = g, item = cl, n_ref = length(ref), n_test = length(mine),
               n_shared = length(intersect(ref, mine)))
  }))
}))
print(par)
stopifnot(par[, all(n_ref == n_test & n_ref == n_shared)])
}

cnt = rbindlist(lapply(c("kept", "removed"), function(mh) {
  x = if (mh == "removed") truth[multiallelic == FALSE] else truth
  x[, .(multiallelic = mh, eligible = .N, informative_negative = sum(informative_negative),
        pos_t25 = sum(pos_t25), pos_t50 = sum(pos_t50), pos_t75 = sum(pos_t75),
        pos_t50_TC_or_T = sum(pos_t50 & type %in% c("TC", "T:C"))), by = .(variant, gate)]
}))
cnt = merge(cnt, truth[, .(multiallelic_eligible = sum(multiallelic), multiallelic_pos_t50 = sum(multiallelic & pos_t50),
                           multiallelic_neg = sum(multiallelic & informative_negative)), by = .(variant, gate)],
            by = c("variant", "gate"))
print(cnt)
fwrite(cnt, file.path(OUT, "truth_counts.tsv"), sep = "\t")
rm(U); invisible(gc())
if (TRUTH_ONLY) quit(save = "no")

# ---- calls (published state: both window filters) -------------------------------------------------
RES = c(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6")
UNIT_RE = c(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
TIER_OUT = c(arm = "dataset", sample = "sample", celltype = "celltype")
COMBO = c("caller", "bam_type", "strand_mode", "snp_subtract")
unreported = function(caller, bam_type) bam_type == "sample" |
  (caller %in% c("scallele", "cellsnp") & bam_type %in% c("recal", "recalq20"))
files = rbindlist(lapply(names(RES), function(tr) {
  rd = file.path(PM, RES[[tr]])
  units = grep(UNIT_RE[[tr]], list.dirs(rd, recursive = FALSE, full.names = FALSE), value = TRUE)
  rbindlist(lapply(units, function(u) {
    f = list.files(file.path(rd, u, "04_filtered"), "_pre100\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
    p = tstrsplit(sub(paste0("^", file.path(rd, u, "04_filtered"), "/"), "", f), "/", fixed = TRUE)
    data.table(tier = tr, unit = u, caller = p[[1]], bam_type = p[[2]], strand_mode = p[[3]],
               snp_subtract = sub("^snp_", "", p[[4]]), file = f)
  }))
}))
files = files[!unreported(caller, bam_type)]
files[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
stopifnot(files[, uniqueN(config), by = tier][, all(V1 == 56 + 2 * any(files$caller == "cellsnp"))])
msg("call files: ", nrow(files), " (", paste(files[, .N, by = tier][, paste(tier, N)], collapse = ", "), ")")
calls = rbindlist(lapply(seq_len(nrow(files)), function(i) {
  if (file.size(files$file[i]) == 0) return(NULL)
  x = fread(files$file[i], select = c("site", "varType"), showProgress = FALSE)[varType %in% c("AG", "TC")]
  unique(x)[, `:=`(tier = files$tier[i], unit = files$unit[i], config = files$config[i])]
}))
msg("A>G / T>C call rows: ", nrow(calls))

# ---- scoring -----------------------------------------------------------------------------------
SCEN = data.table(scenario = c("AG_truth__AG_calls", "AG_truth__AGTC_calls", "AGTC_truth__AGTC_calls", "pooled_truth__AGTC_calls"),
                  variant = c("AG", "AG", "AGTC", "pooled"), call_match = c("AG", "AGTC", "AGTC", "AGTC"))
grid_f1 = function(p, r) fifelse(!is.na(p) & p + r > 0, 2 * p * r / (p + r), NA_real_)
CALLS = list(AG = unique(calls[varType == "AG", .(tier, unit, config, site)]),
             AGTC = unique(calls[, .(tier, unit, config, site)]))
ts = rbindlist(lapply(seq_len(nrow(SCEN)), function(si) rbindlist(lapply(c("kept", "removed"), function(mh)
  rbindlist(lapply(names(GATES), function(g) rbindlist(lapply(names(DS), function(tn) {
    tr = truth[variant == SCEN$variant[si] & gate == g]
    if (mh == "removed") tr = tr[multiallelic == FALSE]
    possite = tr[get(paste0("pos_", tn)) == TRUE, site]; negsite = tr[informative_negative == TRUE, site]
    cl = copy(CALLS[[SCEN$call_match[si]]])
    cl[, `:=`(tp = site %chin% possite, fp = site %chin% negsite)]
    us = cl[, .(n_tp = sum(tp), n_fp = sum(fp)), by = .(tier, unit, config)]
    us = merge(unique(files[, .(tier, unit, config)]), us, by = c("tier", "unit", "config"), all.x = TRUE)
    us[is.na(n_tp), `:=`(n_tp = 0L, n_fp = 0L)]
    us[, `:=`(precision = fifelse(n_tp + n_fp > 0, n_tp / (n_tp + n_fp), NA_real_), recall = n_tp / length(possite))]
    us[, .(precision = mean(precision, na.rm = TRUE), recall = mean(recall), n_tp = mean(n_tp), n_fp = mean(n_fp)),
       by = .(tier, config)][, `:=`(scenario = SCEN$scenario[si], multiallelic = mh, gate = g, t = tn,
                                    n_pos = length(possite), n_neg = length(negsite))]
  }))))))))
ts[is.nan(precision), precision := NA_real_]
ts[, (COMBO) := tstrsplit(config, "__", fixed = TRUE)]
ts[, f1_grid := grid_f1(precision, recall)]
fwrite(ts, file.path(OUT, "tier_scores.tsv"), sep = "\t")

# parity with the published tier scoring
pp = rbindlist(lapply(names(GATES), function(g) rbindlist(lapply(names(DS), function(tn) {
  pub = fread(file.path(V6, "scoring", GATES[[g]], "teratoma_inputs", "invariance", sprintf("d%.2f", DS[[tn]]),
                        "gt_tier_scoring.tsv"))
  pub[, tier := names(TIER_OUT)[match(tier, TIER_OUT)]]
  pub = pub[!unreported(caller, bam_type)]
  me = ts[scenario == "AG_truth__AG_calls" & multiallelic == "kept" & gate == g & t == tn]
  x = merge(me, pub, by = c("tier", COMBO), suffixes = c("", ".pub"))
  data.table(check = "tier_scores", gate = g, item = tn, n_ref = nrow(pub), n_test = nrow(me), n_shared = nrow(x),
             max_abs_diff_precision = max(abs(x$precision - x$precision.pub), na.rm = TRUE),
             max_abs_diff_recall = max(abs(x$recall - x$recall_vs_union)),
             na_mismatch = sum(is.na(x$precision) != is.na(x$precision.pub)))
}))))
print(pp)
fwrite(rbindlist(list(par, pp), fill = TRUE), file.path(OUT, "parity.tsv"), sep = "\t")
stopifnot(pp[, all(n_ref == n_shared & n_test == n_shared & max_abs_diff_precision < 1e-9 &
                   max_abs_diff_recall < 1e-9 & na_mismatch == 0)])

# ---- caller rank (cross-tier rank sum) ---------------------------------------------------------
KEY = c("scenario", "multiallelic", "gate", "t")
ts[, f1_rank := frank(-f1_grid, ties.method = "min"), by = c(KEY, "tier")]
rs = ts[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1_grid, na.rm = TRUE)), by = c(KEY, COMBO)]
rs[, overall_rank := frank(list(rank_sum, -mean_f1), ties.method = "min"), by = KEY]
cb = rs[order(overall_rank), .SD[1], by = c(KEY, "caller")]
cb[, caller_rank := frank(overall_rank, ties.method = "min"), by = KEY]
fwrite(cb, file.path(OUT, "caller_rank.tsv"), sep = "\t")
print(dcast(cb[multiallelic == "kept" & t == "t50"], gate + caller ~ scenario, value.var = "caller_rank"))

sm = ts[, .(recall = median(recall), precision = median(precision, na.rm = TRUE), f1_grid = median(f1_grid, na.rm = TRUE)),
        by = c(KEY, "tier", "caller", "strand_mode")]
fwrite(sm, file.path(OUT, "strand_mode_recall.tsv"), sep = "\t")
print(dcast(sm[multiallelic == "kept" & t == "t50" & gate == "expected" & tier == "arm"],
            caller + strand_mode ~ scenario, value.var = "recall"), digits = 3)
msg("done")
