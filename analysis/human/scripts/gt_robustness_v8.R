#!/usr/bin/env Rscript
# =============================================================================
# gt_robustness_v8.R   (2026-10-13)
#
# One-threshold-at-a-time robustness of the v8 ground truth (gt_eval_v8/README.md), both gates.
# The v6 lattice (gt_robustness_grid_v6.R) ran on v6 positions (strand-agnostic pileups, T>C,
# multiallelic and PGP1 positions included, no nearby-mismatch filter). Here every combination is
# evaluated on the v8 site set:
#   sites       RNA A>G in 10x read orientation (spectrum12 pileups, as gt_truth_agtc_test_v6.R AG)
#   removed     multiallelic (merged NTC third allele >= 3 reads and >= 5% of depth), PGP1 SNV
#               positions, nearby-mismatch fail (>= 2 mismatch types within w1 = 200 on the read
#               strand among the gate's 12-type gate sites without PGP1; gt_truth_v7_variants.R).
#               The cleaning is held fixed (frozen-rule gate sites as neighbours) so only the
#               threshold under examination moves.
#   rules       as gt_truth_agtc_test_v6.R classify(); gate expected: KO cov x NTC edit >= m;
#               gate depth: KO cov >= c. Centre = FROZEN_RULES.tsv, d = 0.5.
#   axes        min_cov 8-12, min_cov_teratoma 8-12, ko_min_expected 1-5 (expected) or
#               ko_min_cov 8-12 (depth), min_edit_reads 2-4, min_ctrl_edit 0.03-0.07,
#               depletion 0.3-0.7, rule strict2 / evaluable / strict3; 1- and 2-step moves
# Parity: centre at d = 0.25 / 0.50 / 0.75 must equal the v8 positives (gt_eval_v8/truth).
#
# Inputs : gt_eval_v6/spectrum12/{positions_noSNP.bed, candidates_all.tsv.gz, pileups/,
#          eval{,_depth}/truth12_sites.tsv.gz}, pgp1_snps/PGP1_hu43860C_hg38_snvs.bed,
#          gt_eval_v8/truth/{gt_v7_truth_sites.tsv.gz, support/nomm_w1_grid.tsv}
# Outputs: gt_eval_v8/robustness/
#          centre_neighbours.tsv   gate x axis x step: n_true, n_neighbour, Jaccard, frac_lost
#          m_curve.tsv             expected gate: sites vs m (0-6), other thresholds at the centre
#          ko_cov_curve.tsv        depth gate: sites vs KO coverage (0-30)
#          cleaning_counts.tsv     AG candidate rows removed by each cleaning step
#          parity.tsv
# Usage  : qsub human/pm/scripts/gt_robustness_v8.sh
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
V6 = file.path(PM, "gt_eval_v6"); SP = file.path(V6, "spectrum12")
V8 = file.path(PM, "gt_eval_v8")
OUT = file.path(V8, "robustness"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

fr = fread(file.path(V6, "FROZEN_RULES.tsv")); P = setNames(fr$value, fr$parameter)
stopifnot(P[["consistency_mode"]] == "evaluable", P[["bq"]] == "30", as.numeric(P[["min_cov_ko"]]) == 0)
KO = c("ADAR1_KO", "ADARB1_KO"); TT = c("T1", "T2", "T3", "T4")
UNITS = c("merged_NTC_ctrl", paste0("merged_", KO), paste0(TT, "_NTC_ctrl"))
COMP = c(A = "T", C = "G", G = "C", T = "A"); BASES = c("A", "C", "G", "T")
ALPHA = 0.05; MIN_THIRD_READS = 3; MIN_THIRD_FRAC = 0.05
GATES = c(expected = "eval", depth = "eval_depth")
KO_AXIS = c(expected = "ko_min_expected", depth = "ko_min_cov")
AXES = list(min_cov = c(8, 9, 10, 11, 12), min_cov_teratoma = c(8, 9, 10, 11, 12),
            ko_min_expected = c(1, 2, 3, 4, 5), ko_min_cov = c(8, 9, 10, 11, 12),
            min_edit_reads = c(2, 3, 4), min_ctrl_edit = c(0.03, 0.04, 0.05, 0.06, 0.07),
            depletion = c(0.3, 0.4, 0.5, 0.6, 0.7), rule = c("strict2", "evaluable", "strict3"))
CENTRE = list(min_cov = as.numeric(P[["min_cov"]]), min_cov_teratoma = as.numeric(P[["min_cov_teratoma"]]),
              ko_min_expected = as.numeric(P[["ko_min_expected"]]), ko_min_cov = 10,
              min_edit_reads = as.numeric(P[["min_edit_reads"]]), min_ctrl_edit = as.numeric(P[["min_ctrl_edit"]]),
              depletion = as.numeric(P[["depletion"]]), rule = "evaluable")

# ---- strand-aware pileups at RNA A>G candidate rows (gt_truth_agtc_test_v6.R) ---------------------
pos = fread(file.path(SP, "positions_noSNP.bed"), header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
pos[, idx := .I]; setkey(pos, chr, pos); N = nrow(pos); msg("positions: ", N)
cand = fread(cmd = paste("gzip -cd", shQuote(file.path(SP, "candidates_all.tsv.gz"))), sep = "\t", header = FALSE,
             select = 1:5, col.names = c("chr", "pos", "strand", "rna_ref", "rna_alt"), showProgress = FALSE)
cand = unique(pos[cand, on = .(chr, pos), nomatch = NULL][rna_ref == "A" & rna_alt == "G", .(chr, pos, idx, strand)])
msg("AG candidate rows: ", nrow(cand))
cols = c("ref_plus", "ref_minus", "A", "C", "G", "T", "a", "c", "g", "t")
plus = cand$strand == "+"
col_of = function(b) ifelse(plus, match(b, cols), match(tolower(COMP[[b]]), cols))
rnd4 = function(x) as.numeric(sprintf("%.4f", x))
S = list(); third = NULL
for (u in UNITS) {
  fs = list.files(file.path(SP, "pileups"), pattern = paste0("^", u, "__.*\\.tsv\\.gz$"), full.names = TRUE)
  d = rbindlist(lapply(fs, function(f) fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", showProgress = FALSE)))
  d = pos[d, on = .(chr, pos), nomatch = NULL]
  m = matrix(0L, N, length(cols), dimnames = list(NULL, cols))
  for (cn in cols) m[d$idx, cn] = as.integer(d[[cn]])
  alt = m[cbind(cand$idx, col_of("G"))]
  refn = ifelse(plus, m[cand$idx, "ref_plus"], m[cand$idx, "ref_minus"])
  cov = refn + alt
  S[[u]] = list(cov = cov, alt = alt, edit = ifelse(cov > 0, rnd4(alt / pmax(cov, 1)), 0))
  if (u == "merged_NTC_ctrl") {
    t1 = m[cbind(cand$idx, col_of("C"))]; t2 = m[cbind(cand$idx, col_of("T"))]
    third = list(mx = pmax(t1, t2), tot = cov + t1 + t2)
  }
  rm(m, d); invisible(gc()); msg(u, ": ", length(fs), " chunks")
}
multi = third$mx >= MIN_THIRD_READS & third$mx >= MIN_THIRD_FRAC * third$tot

# ---- fixed v8 cleaning ---------------------------------------------------------------------------
pgp1 = unique(fread(file.path(PM, "pgp1_snps", "PGP1_hu43860C_hg38_snvs.bed"), header = FALSE, select = c(1, 3),
                    col.names = c("chr", "pos")))
is_pgp1 = paste(cand$chr, cand$pos) %in% paste(pgp1$chr, pgp1$pos)
w1g = fread(file.path(V8, "truth", "support", "nomm_w1_grid.tsv"))[pgp1 == TRUE & chosen == TRUE]
n_types_near = function(q, nb, w) {
  q = copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h = nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL, .(row = i.row, type = x.type)]
  h = rbind(h, q[, .(row, type)])
  r = h[, .(n_types = uniqueN(type)), by = row]
  out = rep(1L, nrow(q)); out[r$row] = r$n_types; out
}
keep = list(); cc = list()
for (g in names(GATES)) {
  t12 = fread(file.path(SP, GATES[[g]], "truth12_sites.tsv.gz"))
  nb = t12[!paste(chr, pos) %in% paste(pgp1$chr, pgp1$pos), .(chr, pos, strand, type)]
  w1 = w1g[gate == g, w1]; stopifnot(length(w1) == 1)
  nomm = n_types_near(cand[, .(chr, pos, strand, type = "AG")], nb, w1) < 2
  keep[[g]] = !multi & !is_pgp1 & nomm
  cc[[g]] = data.table(gate = g, w1 = w1, n_ag_rows = nrow(cand), multiallelic = sum(multi), pgp1 = sum(is_pgp1),
                       nomm_fail = sum(!nomm), kept = sum(keep[[g]]))
  msg(g, ": w1 = ", w1, "; kept ", sum(keep[[g]]), " of ", nrow(cand))
}
fwrite(rbindlist(cc), file.path(OUT, "cleaning_counts.tsv"), sep = "\t")

# ---- rules ---------------------------------------------------------------------------------------
cp0 = function(n) { o = rep(1, length(n)); ok = n > 0; o[ok] = qbeta(1 - ALPHA, 1, n[ok]); o }
UB = lapply(KO, function(a) cp0(S[[paste0("merged_", a)]]$cov))
truth_of = function(g, p) {
  ctrl = S[["merged_NTC_ctrl"]]
  cov_arm = lapply(KO, function(a) { k = S[[paste0("merged_", a)]]
    if (g == "expected") ctrl$cov >= p$min_cov & k$cov * ctrl$edit >= p$ko_min_expected & (p$ko_min_expected > 0 | k$cov > 0)
    else ctrl$cov >= p$min_cov & k$cov >= p$ko_min_cov & k$cov > 0 })
  ctrl_ed = ctrl$alt >= p$min_edit_reads & ctrl$edit >= p$min_ctrl_edit
  n_eval = Reduce(`+`, lapply(TT, function(t) S[[paste0(t, "_NTC_ctrl")]]$cov >= p$min_cov_teratoma))
  n_pos = Reduce(`+`, lapply(TT, function(t) { x = S[[paste0(t, "_NTC_ctrl")]]
    x$cov >= p$min_cov_teratoma & x$alt >= p$min_edit_reads & x$edit >= p$min_ctrl_edit }))
  cons = switch(p$rule, strict2 = n_pos >= 2, strict3 = n_pos >= 3, evaluable = n_pos >= 2 & n_pos >= pmax(2, n_eval - 1))
  gate = Reduce(`|`, cov_arm) & ctrl_ed & cons & keep[[g]]
  target = (1 - p$depletion) * ctrl$edit
  dep = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & S[[paste0("merged_", KO[i])]]$edit <= target))
  pw = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & UB[[i]] < target))
  list(gated = which(gate), true = which(gate & dep), powered = which(gate & dep & pw))
}

# ---- parity: centre at the three t equals the v8 positives -------------------------------------
v8 = fread(file.path(V8, "truth", "gt_v7_truth_sites.tsv.gz"))
key = paste(cand$chr, cand$pos, cand$strand)
par = rbindlist(lapply(names(GATES), function(g) rbindlist(lapply(c(t25 = 0.25, t50 = 0.50, t75 = 0.75), function(d) {
  p = modifyList(CENTRE, list(depletion = d)); tr = key[truth_of(g, p)$true]
  tn = c(`0.25` = "t25", `0.5` = "t50", `0.75` = "t75")[[as.character(d)]]
  ref = v8[gate == g & get(paste0("class_", tn)) == "positive", paste(chr, pos, read_strand)]
  data.table(gate = g, t = tn, n_lattice = length(tr), n_v8 = length(ref), n_shared = length(intersect(tr, ref)))
}))))
par[, identical := n_lattice == n_v8 & n_shared == n_v8]
print(par); fwrite(par, file.path(OUT, "parity.tsv"), sep = "\t")
stopifnot(par$identical)
msg("PARITY OK: centre reproduces the v8 positives at t25 / t50 / t75, both gates")

# ---- centre neighbours ---------------------------------------------------------------------------
nb_all = list()
for (g in names(GATES)) {
  axes = setdiff(names(AXES), setdiff(KO_AXIS, KO_AXIS[[g]]))
  c0 = truth_of(g, CENTRE)$true
  nb_all[[length(nb_all) + 1]] = data.table(gate = g, axis = axes, step = 0L,
    neighbour_value = vapply(axes, function(a) as.character(CENTRE[[a]]), ""), n_true = length(c0),
    n_neighbour = length(c0), jaccard = 1, frac_lost = 0)
  for (a in axes) {
    i0 = match(as.character(CENTRE[[a]]), as.character(AXES[[a]]))
    for (s in c(-2, -1, 1, 2)) {
      j = i0 + s; if (j < 1 || j > length(AXES[[a]])) next
      p = CENTRE; p[[a]] = AXES[[a]][[j]]
      tr = truth_of(g, p)$true; inter = length(intersect(c0, tr))
      nb_all[[length(nb_all) + 1]] = data.table(gate = g, axis = a, step = as.integer(s),
        neighbour_value = as.character(AXES[[a]][[j]]), n_true = length(c0), n_neighbour = length(tr),
        jaccard = inter / (length(c0) + length(tr) - inter), frac_lost = 1 - inter / length(c0))
    }
  }
  msg(g, ": centre ", length(c0), " sites")
}
nb_all = rbindlist(nb_all); setorder(nb_all, gate, axis, step)
fwrite(nb_all, file.path(OUT, "centre_neighbours.tsv"), sep = "\t")

# ---- informative-KO curves -------------------------------------------------------------------------
curve = function(g, axis, vals) rbindlist(lapply(vals, function(v) {
  p = CENTRE; p[[axis]] = v; r = truth_of(g, p)
  n_t = vapply(c(0, 0.25, 0.75), function(d) length(truth_of(g, modifyList(p, list(depletion = d)))$true), 1L)
  data.table(value = v, n_gated = length(r$gated), n_true = length(r$true), n_true_powered = length(r$powered),
             n_true_t0 = n_t[1], n_true_t25 = n_t[2], n_true_t75 = n_t[3])
}))
mc = curve("expected", "ko_min_expected", sort(unique(round(c(seq(0, 6, 0.1), -log(0.05), -log(0.01)), 4))))
setnames(mc, "value", "ko_min_expected"); fwrite(mc, file.path(OUT, "m_curve.tsv"), sep = "\t")
kc = curve("depth", "ko_min_cov", 0:30)
setnames(kc, "value", "ko_min_cov"); fwrite(kc, file.path(OUT, "ko_cov_curve.tsv"), sep = "\t")
msg("m curve: TRUE at m = 3: ", mc[abs(ko_min_expected - 3) < 1e-9, n_true],
    "; KO cov curve: TRUE at 10: ", kc[ko_min_cov == 10, n_true])

# ---- optional: full truth tables for the +-1-step neighbours (rank sensitivity, user 2026-10-14) -------
# GT_ROBUST_TRUTH_OUT=<dir> writes <dir>/<axis>_<value>/gt_v7_truth_sites.tsv.gz in the v8 truth format
# (classes as gt_truth_agtc_test_v6.R classify(): informative negative = both KO arms covered with CP lower
# bound > 75% of NTC editing). Depletion is the t25 / t50 / t75 axis of every truth, so it is not swept here.
# A gate-specific KO axis keeps the other gate at the centre. Parity: the centre table equals the v8 truth.
TRUTH_OUT = Sys.getenv("GT_ROBUST_TRUTH_OUT", "")
if (nzchar(TRUTH_OUT)) {
  DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75); NEG_MAX_DEP = 0.25
  cp_lower = function(k, n) { o = rep(0, length(k)); ok = n > 0 & k > 0; o[ok] = qbeta(ALPHA, k[ok], n[ok] - k[ok] + 1); o }
  classes_of = function(g, p) {
    ctrl = S[["merged_NTC_ctrl"]]; gated = truth_of(g, p)$gated
    cov_arm = lapply(KO, function(a) { k = S[[paste0("merged_", a)]]
      if (g == "expected") ctrl$cov >= p$min_cov & k$cov * ctrl$edit >= p$ko_min_expected & (p$ko_min_expected > 0 | k$cov > 0)
      else ctrl$cov >= p$min_cov & k$cov >= p$ko_min_cov & k$cov > 0 })
    nd = Reduce(`&`, lapply(seq_along(KO), function(i) { k = S[[paste0("merged_", KO[i])]]
      cov_arm[[i]] & cp_lower(k$alt, k$cov) > (1 - NEG_MAX_DEP) * ctrl$edit }))
    o = data.table(gate = g, site = paste0(cand$chr[gated], ":", cand$pos[gated]), chr = cand$chr[gated],
                   pos = cand$pos[gated], genomic_ref = fifelse(cand$strand[gated] == "+", "A", "T"),
                   a_to_i_strand = cand$strand[gated], read_type = "AG", read_strand = cand$strand[gated])
    for (tn in names(DS)) {
      pt = gated %in% truth_of(g, modifyList(p, list(depletion = DS[[tn]])))$true
      stopifnot(!any(pt & nd[gated]))
      set(o, j = paste0("class_", tn), value = fifelse(pt, "positive", fifelse(nd[gated], "informative_negative", "unresolved")))
    }
    o
  }
  write_truth = function(lab, axis = NULL, value = NULL) {
    tr = rbindlist(lapply(names(GATES), function(g) {
      p = CENTRE
      if (!is.null(axis) && (!axis %in% KO_AXIS || KO_AXIS[[g]] == axis)) p[[axis]] = value
      classes_of(g, p)
    }))
    setorder(tr, gate, chr, pos)
    d = file.path(TRUTH_OUT, lab); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    fwrite(tr, file.path(d, "gt_v7_truth_sites.tsv.gz"), sep = "\t")
    tr[, .(truth = lab, n_pos_t25 = sum(class_t25 == "positive"), n_pos_t50 = sum(class_t50 == "positive"),
           n_pos_t75 = sum(class_t75 == "positive"), n_neg = sum(class_t50 == "informative_negative"), n_sites = .N), by = gate]
  }
  cnt = list(write_truth("centre"))
  ref = fread(file.path(V8, "truth", "gt_v7_truth_sites.tsv.gz"))
  cen = fread(file.path(TRUTH_OUT, "centre", "gt_v7_truth_sites.tsv.gz"))
  k = c("gate", "site", "read_strand", paste0("class_", names(DS)))
  stopifnot(nrow(cen) == nrow(ref), fsetequal(cen[, ..k], ref[, ..k]))
  msg("TRUTH PARITY OK: centre table equals the v8 truth (", nrow(ref), " rows)")
  for (a in setdiff(names(AXES), "depletion")) {
    i0 = match(as.character(CENTRE[[a]]), as.character(AXES[[a]]))
    for (s in c(-1, 1)) {
      j = i0 + s; if (j < 1 || j > length(AXES[[a]])) next
      cnt[[length(cnt) + 1]] = write_truth(paste0(a, "_", AXES[[a]][[j]]), a, AXES[[a]][[j]])
    }
  }
  cnt = rbindlist(cnt); print(cnt); fwrite(cnt, file.path(TRUTH_OUT, "truth_counts.tsv"), sep = "\t")
}
msg("done -> ", OUT)
