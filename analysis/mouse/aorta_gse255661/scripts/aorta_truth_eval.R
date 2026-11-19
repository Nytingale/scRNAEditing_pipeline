#!/usr/bin/env Rscript
# =============================================================================
# aorta_truth_eval.R
#
# GSE255661 aorta ground truth (3 Adar1KD vs 3 Control mice, Myh11+ cells, 10 weeks): the teratoma v8
# truth rules (human/pm/gt_eval_v6/FROZEN_RULES.tsv, gt_spectrum12_eval_v6.R, gt_truth_agtc_test_v6.R,
# gt_truth_v7.R, gt_truth_v7_variants.R variant nomm) with one KO arm and the 3 control mice as replicates.
#   12-type gate (per type R>Y, read orientation = RNA strand, 10x 3' forward):
#     covered      expected gate: CTRL cov >= 10 & KD cov x CTRL edit >= 3;  depth gate: CTRL cov >= 10 & KD cov >= 10
#     ctrl_edited  merged CTRL alt >= 3 & edit >= 0.05
#     consistent   per control mouse: evaluable = cov >= 10, positive = evaluable & alt >= 3 & edit >= 0.05;
#                  n_pos >= 2 & n_pos >= max(2, n_eval - 1)
#     d_max        1 - KD edit / CTRL edit;  positive(t) = gate & d_max >= t  (t = 25 / 50 / 75%)
#     informative negative  gate & one-sided Clopper-Pearson lower bound (alpha 0.05) on KD edit > 0.75 x CTRL edit,
#                  and NOT in REDIportal. Unlike the teratoma (not depleted in ADAR1 KO AND ADARB1 KO), the aorta has
#                  only an Adar1 arm, so "not depleted" also holds at ADAR2 sites (Copa I164V: 50% CTRL, 55% KD);
#                  a call there is correct editing. Negatives in REDIportal are therefore unscored
#                  (class unscored_negative_rediportal; user 2026-10-26).
#   A-to-G truth   type AG gate sites; multiallelic removed (third allele >= 3 reads and >= 5% of depth in merged
#                  CTRL on the read strand); removal of sites near other mismatch types: >= 2 mismatch types within
#                  +- w1 nt on the read strand among all 12-type gate sites of the gate, w1 searched 20-200 by 20 and
#                  chosen by the highest A-to-G share of the 12-type t50 positives that remain (first on ties).
#   SNP modes      excluded = MGP (all strains) SNP positions removed before anything else (also as neighbours);
#                  included = kept. Both are written.
# Input : truth/spectrum12/ (aorta_truth_spectrum12.sh), truth/design.tsv
# Output: truth/eval/<snp_mode>/{truth12_sites.tsv.gz, spectrum12_by_d.tsv, nomm_w1_grid.tsv,
#           gt_aorta_truth_sites.tsv.gz, truth_counts.tsv, positives_by_kd_mouse.tsv,
#           unit_cov.tsv.gz (unit CTRL + control mice, for gt_score_v7.R GT_UNIT_COV)}, truth/eval/summary.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
A = file.path(Sys.getenv("SCRNA_DIR", "."), "mouse/aorta_gse255661")
SP = Sys.getenv("AORTA_SP12_DIR", file.path(A, "truth/spectrum12"))
OUT0 = Sys.getenv("AORTA_EVAL_DIR", file.path(A, "truth/eval")); dir.create(OUT0, showWarnings = FALSE)
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

MIN_COV = 10; MIN_COV_KO_DEPTH = 10; M_EXP = 3; K = 3; E = 0.05; TCOV = 10
NEG_MAX_DEP = 0.25; CI_ALPHA = 0.05; MIN_THIRD_READS = 3; MIN_THIRD_FRAC = 0.05
DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
W_ON = seq(20, 200, 20)
GATES = c("expected", "depth")

design = fread(file.path(A, "truth/design.tsv"))
CM = design[arm == "CTRL", sample]; KM = design[arm == "KD", sample]
stopifnot(length(CM) == 3, length(KM) == 3)
UNITS = c("merged_CTRL", "merged_KD", CM, KM)
COMP = c(A = "T", C = "G", G = "C", T = "A")
BASES = c("A", "C", "G", "T")
TYPES = CJ(R = BASES, Y = BASES)[R != Y][, type := paste0(R, Y)]

# ---- positions and raw counts --------------------------------------------------------------------
rd_bed = function(f) fread(f, header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
pos = rd_bed(file.path(SP, "positions_all.bed"))[, idx := .I]
setkey(pos, chr, pos)
pos[rd_bed(file.path(SP, "positions_mgpSNP.bed")), snp := TRUE, on = .(chr, pos)][is.na(snp), snp := FALSE]
pos[rd_bed(file.path(SP, "positions_rediportal.bed")), redi := TRUE, on = .(chr, pos)][is.na(redi), redi := FALSE]
N = nrow(pos); msg("positions: ", N, " (MGP SNP ", sum(pos$snp), ", REDIportal ", sum(pos$redi), ")")

cols = c("ref_plus", "ref_minus", "A", "C", "G", "T", "a", "c", "g", "t")
nreg = nrow(fread(file.path(SP, "regions.tsv"), header = FALSE))
read_unit = function(u) {
  fs = list.files(file.path(SP, "pileups"), pattern = paste0("^", u, "__.*\\.tsv\\.gz$"), full.names = TRUE)
  stopifnot(length(fs) == nreg)
  d = rbindlist(lapply(fs, function(f) fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", showProgress = FALSE)))
  d = pos[d, on = .(chr, pos), nomatch = NULL]
  m = matrix(0L, N, length(cols), dimnames = list(NULL, cols))
  for (cn in cols) m[d$idx, cn] = as.integer(d[[cn]])
  ref = rep(NA_character_, N); ref[d$idx] = toupper(d$ref)
  msg(u, ": ", nrow(d), " covered positions from ", length(fs), " chunks")
  list(m = m, ref = ref)
}
U = setNames(lapply(UNITS, read_unit), UNITS)
REF = U[["merged_CTRL"]]$ref
for (u in UNITS) REF[is.na(REF)] = U[[u]]$ref[is.na(REF)]

cand = fread(cmd = paste("gzip -cd", shQuote(file.path(SP, "candidates_all.tsv.gz"))), sep = "\t", header = FALSE,
             select = 1:5, col.names = c("chr", "pos", "strand", "rna_ref", "rna_alt"), showProgress = FALSE)
cand = pos[cand, on = .(chr, pos), nomatch = NULL][, type := paste0(rna_ref, rna_alt)]
msg("candidate (position, strand, alt) rows: ", nrow(cand))

rnd4 = function(x) as.numeric(sprintf("%.4f", x))
counts_for = function(u, idx, strand, Y) {
  m = U[[u]]$m; plus = strand == "+"
  alt = ifelse(plus, m[cbind(idx, match(Y, cols))], m[cbind(idx, match(tolower(COMP[[Y]]), cols))])
  cov = ifelse(plus, m[idx, "ref_plus"], m[idx, "ref_minus"]) + alt
  list(cov = cov, alt = alt, edit = ifelse(cov > 0, rnd4(alt / pmax(cov, 1)), 0))
}
third_max = function(idx, strand, R, Y) {
  m = U[["merged_CTRL"]]$m; plus = strand == "+"; th = setdiff(BASES, c(R, Y))
  f = function(b) ifelse(plus, m[cbind(idx, match(b, cols))], m[cbind(idx, match(tolower(COMP[[b]]), cols))])
  list(t1 = f(th[1]), t2 = f(th[2]))
}
cp_lower = function(k, n, a = CI_ALPHA) { o = rep(0, length(k)); ok = n > 0 & k > 0; o[ok] = qbeta(a, k[ok], n[ok] - k[ok] + 1); o }

# ---- 12-type gate sites per gate (all positions; SNP flag carried) ---------------------------------
sites12 = list(); by_d = list()
for (ti in seq_len(nrow(TYPES))) {
  R = TYPES$R[ti]; Y = TYPES$Y[ti]; ty = TYPES$type[ti]
  cd = cand[type == ty]
  if (!nrow(cd)) next
  gref = ifelse(cd$strand == "+", R, COMP[[R]])
  stopifnot(all(REF[cd$idx] == gref, na.rm = TRUE))
  S = lapply(setNames(UNITS, UNITS), function(u) counts_for(u, cd$idx, cd$strand, Y))
  ctrl = S$merged_CTRL; kd = S$merged_KD
  ctrl_ed = ctrl$alt >= K & ctrl$edit >= E
  n_eval = Reduce(`+`, lapply(CM, function(t) S[[t]]$cov >= TCOV))
  n_pos = Reduce(`+`, lapply(CM, function(t) { x = S[[t]]; x$cov >= TCOV & x$alt >= K & x$edit >= E }))
  cons = n_pos >= 2 & n_pos >= pmax(2, n_eval - 1)
  th = third_max(cd$idx, cd$strand, R, Y); mx = pmax(th$t1, th$t2)
  multi = mx >= MIN_THIRD_READS & mx >= MIN_THIRD_FRAC * (ctrl$cov + th$t1 + th$t2)
  for (g in GATES) {
    covd = if (g == "expected") ctrl$cov >= MIN_COV & kd$cov * ctrl$edit >= M_EXP else ctrl$cov >= MIN_COV & kd$cov >= MIN_COV_KO_DEPTH
    gate = covd & ctrl_ed & cons
    neg = gate & cp_lower(kd$alt, kd$cov) > (1 - NEG_MAX_DEP) * ctrl$edit
    w = which(gate)
    if (!length(w)) next
    dt = data.table(gate = g, type = ty, chr = cd$chr[w], pos = cd$pos[w], strand = cd$strand[w],
                    snp = cd$snp[w], redi = cd$redi[w], multiallelic = multi[w],
                    ctrl_cov = ctrl$cov[w], ctrl_alt = ctrl$alt[w], ctrl_edit = ctrl$edit[w],
                    kd_cov = kd$cov[w], kd_alt = kd$alt[w], kd_edit = kd$edit[w],
                    n_ctrl_eval = n_eval[w], n_ctrl_pos = n_pos[w],
                    d_max = 1 - kd$edit[w] / ctrl$edit[w], informative_negative = neg[w])
    for (j in seq_along(KM)) {   # per KD mouse, for the consistency check only (not a rule)
      x = S[[KM[j]]]
      set(dt, j = paste0("kd", j, "_cov"), value = x$cov[w]); set(dt, j = paste0("kd", j, "_alt"), value = x$alt[w])
    }
    for (j in seq_along(CM)) {   # per control mouse: sample-tier recall denominators (unit_cov)
      x = S[[CM[j]]]
      set(dt, j = paste0("c", j, "_cov"), value = x$cov[w]); set(dt, j = paste0("c", j, "_alt"), value = x$alt[w])
    }
    sites12[[paste(g, ty)]] = dt
  }
  msg(ty, ": candidates ", nrow(cd))
}
s12 = rbindlist(sites12)
rm(U); invisible(gc())

# number of mismatch types among `nb` within +-w on the query's strand (downstream_filter filReliable rule)
n_types_near = function(q, nb, w) {
  q = copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h = nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL, .(row = i.row, type = x.type)]
  h = rbind(h, q[, .(row, type)])
  h[, .(n_types = uniqueN(type)), by = row][order(row), n_types]
}

summ = list()
for (sm in c("excluded", "included")) {
  OUT = file.path(OUT0, sm); dir.create(OUT, showWarnings = FALSE)
  t12 = if (sm == "excluded") s12[snp == FALSE] else copy(s12)
  grid = list(); tr = list(); cnt = list()
  for (g in GATES) {
    u = t12[gate == g]
    nb = u[, .(chr, pos, strand, type)]
    p50 = u[d_max >= DS[["t50"]], .(chr, pos, strand, type)]
    gr = rbindlist(lapply(W_ON, function(w) {
      keep = n_types_near(p50, nb, w) < 2
      data.table(gate = g, w1 = w, n_kept = sum(keep), n_ag = sum(keep & p50$type == "AG"))
    }))
    gr[, frac_ag := n_ag / pmax(n_kept, 1)]
    sc = if (uniqueN(gr$frac_ag) == 1) rep(0, nrow(gr)) else as.numeric(scale(gr$frac_ag))
    w1b = gr$w1[which.max(sc)]; gr[, chosen := w1 == w1b]
    grid[[g]] = gr
    t12[gate == g, nomm_pass := n_types_near(u[, .(chr, pos, strand, type)], nb, w1b) < 2]
    t12[gate == g, w1_chosen := w1b]

    ag = t12[gate == g & type == "AG"]
    x = ag[multiallelic == FALSE & nomm_pass == TRUE]
    x[, `:=`(site = paste0(chr, ":", pos), genomic_ref = fifelse(strand == "+", "A", "T"), a_to_i_strand = strand,
             read_type = type, read_strand = strand)]
    # Negatives in REDIportal are left unscored: with only an Adar1 arm, "not depleted" includes ADAR2 sites
    # (e.g. Copa I164V), where a call is correct editing, not a false positive.
    for (tn in names(DS)) set(x, j = paste0("class_", tn),
      value = fifelse(x$d_max >= DS[[tn]], "positive",
              fifelse(x$informative_negative & x$redi, "unscored_negative_rediportal",
              fifelse(x$informative_negative, "informative_negative", "unresolved"))))
    stopifnot(!any(x$informative_negative & x$d_max >= DS[["t25"]]))
    tr[[g]] = x
    for (tn in names(DS)) {
      step = list(gate_AG = ag, multiallelic_removed = ag[multiallelic == FALSE], final = x)
      cnt[[paste(g, tn)]] = rbindlist(lapply(names(step), function(st) { y = step[[st]]
        data.table(snp_mode = sm, gate = g, t = tn, stage = st, w1 = w1b, eligible = nrow(y),
                   positive = sum(y$d_max >= DS[[tn]]), informative_negative = sum(y$informative_negative & !y$redi),
                   negative_in_rediportal_unscored = sum(y$informative_negative & y$redi),
                   positive_in_rediportal = sum(y$redi & y$d_max >= DS[[tn]])) }))
    }
  }
  tr = rbindlist(tr); cnt = rbindlist(cnt)
  COLS = c("gate", "site", "chr", "pos", "genomic_ref", "a_to_i_strand", "read_type", "read_strand", paste0("class_", names(DS)))
  setorder(tr, gate, chr, pos)
  fwrite(tr[, ..COLS], file.path(OUT, "gt_aorta_truth_sites.tsv.gz"), sep = "\t")
  fwrite(t12, file.path(OUT, "truth12_sites.tsv.gz"), sep = "\t")
  # unit coverage at the truth positions (gt_score_v7.R GT_UNIT_COV format): unit CTRL = merged arm, then the mice
  us = unique(tr[, c("site", "ctrl_cov", "ctrl_alt", paste0("c", seq_along(CM), c("_cov")), paste0("c", seq_along(CM), "_alt")), with = FALSE])
  stopifnot(!anyDuplicated(us$site))
  uc = rbindlist(c(list(us[, .(unit = "CTRL", site, coverage = ctrl_cov, variant_cov = ctrl_alt)]),
    lapply(seq_along(CM), function(j) us[, .(unit = CM[j], site, coverage = get(paste0("c", j, "_cov")),
                                              variant_cov = get(paste0("c", j, "_alt")))])))
  fwrite(uc, file.path(OUT, "unit_cov.tsv.gz"), sep = "\t")
  fwrite(rbindlist(grid), file.path(OUT, "nomm_w1_grid.tsv"), sep = "\t")
  fwrite(cnt, file.path(OUT, "truth_counts.tsv"), sep = "\t")
  # 12-type positives by type (before / after removal of sites near other mismatch types)
  sp = rbindlist(lapply(names(DS), function(tn) rbind(
    t12[d_max >= DS[[tn]], .(N = .N), by = .(gate, type)][, `:=`(t = tn, cleaning = "none")],
    t12[d_max >= DS[[tn]] & nomm_pass == TRUE, .(N = .N), by = .(gate, type)][, `:=`(t = tn, cleaning = "nomm")])))
  sp[, pct := 100 * N / sum(N), by = .(gate, t, cleaning)]
  fwrite(sp, file.path(OUT, "spectrum12_by_t.tsv"), sep = "\t")
  # depletion in each KD mouse at the t50 positives (KD edit <= 0.5 x CTRL edit, mice with >= 10 reads)
  pk = rbindlist(lapply(seq_along(KM), function(j) {
    p = tr[class_t50 == "positive"]
    p[, `:=`(cv = get(paste0("kd", j, "_cov")), al = get(paste0("kd", j, "_alt")))]
    p[, .(kd_mouse = KM[j], n_pos = .N, n_cov10 = sum(cv >= 10),
          n_depleted = sum(cv >= 10 & al / pmax(cv, 1) <= 0.5 * ctrl_edit)), by = gate]
  }))
  fwrite(pk, file.path(OUT, "positives_by_kd_mouse.tsv"), sep = "\t")
  summ[[sm]] = cnt
  msg("SNP ", sm, ":"); print(cnt[stage == "final"]); print(pk)
}
fwrite(rbindlist(summ), file.path(OUT0, "summary.tsv"), sep = "\t")
msg("done")
