#!/usr/bin/env Rscript
# =============================================================================
# gt_robustness_v9.R   (2026-11-16; v9 port of gt_robustness_v8.R for Response Figure 1)
#
# One-threshold-at-a-time robustness of the v9 teratoma ground truth (gt_truth_v9.R), both gates, on the
# seed-8787 spectrum12 pileups. As v8 except for the v9 site definition:
#   sites       read-orientation A>G and T>C candidate rows (gt_truth_agtc_test_v6.R AGTC); at each parameter
#               setting, where both orientations of a position are eligible, the one with more merged-NTC reads
#               is kept (ties: A>G), then
#   removed     multiallelic (merged NTC third allele >= 3 reads and >= 5% of depth, on the kept orientation),
#               PGP1 SNV positions, nearby-mismatch fail (>= 2 mismatch types within w1 on the read strand among
#               the gate's 12-type gate sites without PGP1; w1 from gt_eval_v9/build/variants/nomm_w1_grid.tsv)
#   rules       gate expected: KO cov x NTC edit >= m; gate depth: KO cov >= c. Centre = FROZEN_RULES.tsv.
#   axes        as v8 (1- and 2-step moves); curves: m 0-6 (expected), KO coverage 0-30 (depth)
# Parity: centre at d = 0.25 / 0.50 / 0.75 must equal the v9 positives (gt_eval_v9/truth/gt_v9_truth_sites.tsv.gz,
# chr:pos and read strand).
# Inputs : gt_eval_v6/spectrum12/{positions_noSNP.bed, candidates_all.tsv.gz, pileups/},
#          gt_eval_v9/build/variants/{truth12_site_flags.tsv.gz, nomm_w1_grid.tsv}, pgp1_snps/PGP1_hu43860C_hg38_snvs.bed
# Outputs: gt_eval_v9/robustness/{centre_neighbours, m_curve, ko_cov_curve, cleaning_counts, parity}.tsv
#          (columns as v8; read by GTDEF_prep_source_data.R with GT_V9=1)
# Usage  : qsub human/pm/scripts/gt_robustness_v9.sh
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
V6 = file.path(PM, "gt_eval_v6"); SP = file.path(V6, "spectrum12")
V9 = file.path(PM, "gt_eval_v9")
OUT = file.path(V9, "robustness"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

fr = fread(file.path(V6, "FROZEN_RULES.tsv")); P = setNames(fr$value, fr$parameter)
stopifnot(P[["consistency_mode"]] == "evaluable", P[["bq"]] == "30", as.numeric(P[["min_cov_ko"]]) == 0)
KO = c("ADAR1_KO", "ADARB1_KO"); TT = c("T1", "T2", "T3", "T4")
UNITS = c("merged_NTC_ctrl", paste0("merged_", KO), paste0(TT, "_NTC_ctrl"))
COMP = c(A = "T", C = "G", G = "C", T = "A"); BASES = c("A", "C", "G", "T")
ALPHA = 0.05; MIN_THIRD_READS = 3; MIN_THIRD_FRAC = 0.05
GATES = c("expected", "depth")
KO_AXIS = c(expected = "ko_min_expected", depth = "ko_min_cov")
AXES = list(min_cov = c(8, 9, 10, 11, 12), min_cov_teratoma = c(8, 9, 10, 11, 12),
            ko_min_expected = c(1, 2, 3, 4, 5), ko_min_cov = c(8, 9, 10, 11, 12),
            min_edit_reads = c(2, 3, 4), min_ctrl_edit = c(0.03, 0.04, 0.05, 0.06, 0.07),
            depletion = c(0.3, 0.4, 0.5, 0.6, 0.7), rule = c("strict2", "evaluable", "strict3"))
CENTRE = list(min_cov = as.numeric(P[["min_cov"]]), min_cov_teratoma = as.numeric(P[["min_cov_teratoma"]]),
              ko_min_expected = as.numeric(P[["ko_min_expected"]]), ko_min_cov = 10,
              min_edit_reads = as.numeric(P[["min_edit_reads"]]), min_ctrl_edit = as.numeric(P[["min_ctrl_edit"]]),
              depletion = as.numeric(P[["depletion"]]), rule = "evaluable")

# ---- strand-aware pileups at read-orientation A>G and T>C candidate rows ------------------------
pos = fread(file.path(SP, "positions_noSNP.bed"), header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
pos[, idx := .I]; setkey(pos, chr, pos); N = nrow(pos); msg("positions: ", N)
cand = fread(cmd = paste("gzip -cd", shQuote(file.path(SP, "candidates_all.tsv.gz"))), sep = "\t", header = FALSE,
             select = 1:5, col.names = c("chr", "pos", "strand", "rna_ref", "rna_alt"), showProgress = FALSE)
cand = pos[cand, on = .(chr, pos), nomatch = NULL]
cand[, type := paste0(rna_ref, rna_alt)]
cand = unique(cand[type %in% c("AG", "TC"), .(chr, pos, idx, strand, type)])
setorder(cand, chr, pos, type, strand)
msg("AG / TC candidate rows: ", paste(cand[, .N, by = type][, paste(type, N)], collapse = ", "))
cols = c("ref_plus", "ref_minus", "A", "C", "G", "T", "a", "c", "g", "t")
plus = cand$strand == "+"
R_B = substr(cand$type, 1, 1); Y_B = substr(cand$type, 2, 2)
col_of = function(b) ifelse(plus, match(b, cols), match(tolower(COMP[b]), cols))
rnd4 = function(x) as.numeric(sprintf("%.4f", x))
th1 = ifelse(cand$type == "AG", "C", "A"); th2 = ifelse(cand$type == "AG", "T", "G")
S = list(); third = NULL
for (u in UNITS) {
  fs = list.files(file.path(SP, "pileups"), pattern = paste0("^", u, "__.*\\.tsv\\.gz$"), full.names = TRUE)
  d = rbindlist(lapply(fs, function(f) fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", showProgress = FALSE)))
  d = pos[d, on = .(chr, pos), nomatch = NULL]
  m = matrix(0L, N, length(cols), dimnames = list(NULL, cols))
  for (cn in cols) m[d$idx, cn] = as.integer(d[[cn]])
  alt = m[cbind(cand$idx, col_of(Y_B))]
  refn = ifelse(plus, m[cand$idx, "ref_plus"], m[cand$idx, "ref_minus"])
  cov = refn + alt
  S[[u]] = list(cov = cov, alt = alt, edit = ifelse(cov > 0, rnd4(alt / pmax(cov, 1)), 0))
  if (u == "merged_NTC_ctrl") {
    t1 = m[cbind(cand$idx, col_of(th1))]; t2 = m[cbind(cand$idx, col_of(th2))]
    third = list(mx = pmax(t1, t2), tot = cov + t1 + t2)
  }
  rm(m, d); invisible(gc()); msg(u, ": ", length(fs), " chunks")
}
multi = third$mx >= MIN_THIRD_READS & third$mx >= MIN_THIRD_FRAC * third$tot
site_key = paste(cand$chr, cand$pos)

# ---- fixed v9 cleaning (applied after the orientation choice) -------------------------------------
pgp1 = unique(fread(file.path(PM, "pgp1_snps", "PGP1_hu43860C_hg38_snvs.bed"), header = FALSE, select = c(1, 3),
                    col.names = c("chr", "pos")))
is_pgp1 = site_key %in% paste(pgp1$chr, pgp1$pos)
w1g = fread(file.path(V9, "build", "variants", "nomm_w1_grid.tsv"))[pgp1 == TRUE & chosen == TRUE]
t12 = fread(file.path(V9, "build", "variants", "truth12_site_flags.tsv.gz"))[pgp1 == FALSE]
n_types_near = function(q, nb, w) {
  q = copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h = nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL, .(row = i.row, type = x.type)]
  h = rbind(h, q[, .(row, type)])
  r = h[, .(n_types = uniqueN(type)), by = row]
  out = rep(1L, nrow(q)); out[r$row] = r$n_types; out
}
clean = list(); cc = list()
for (g in GATES) {
  nb = t12[gate == g, .(chr, pos, strand, type)]
  w1 = w1g[gate == g, w1]; stopifnot(length(w1) == 1)
  nomm = n_types_near(cand[, .(chr, pos, strand, type)], nb, w1) < 2
  clean[[g]] = !multi & !is_pgp1 & nomm
  cc[[g]] = data.table(gate = g, w1 = w1, n_agtc_rows = nrow(cand), multiallelic = sum(multi), pgp1 = sum(is_pgp1),
                       nomm_fail = sum(!nomm), kept = sum(clean[[g]]))
  msg(g, ": w1 = ", w1, "; clean rows ", sum(clean[[g]]), " of ", nrow(cand))
}
fwrite(rbindlist(cc), file.path(OUT, "cleaning_counts.tsv"), sep = "\t")

# ---- rules ---------------------------------------------------------------------------------------
cp0 = function(n) { o = rep(1, length(n)); ok = n > 0; o[ok] = qbeta(1 - ALPHA, 1, n[ok]); o }
UB = lapply(KO, function(a) cp0(S[[paste0("merged_", a)]]$cov))
ntc_cov = S[["merged_NTC_ctrl"]]$cov
# per position, the eligible orientation with more merged-NTC reads (ties: A>G); returns row indices
pick = function(elig) {
  i = which(elig)
  if (!length(i)) return(integer(0))
  o = data.table(i = i, k = site_key[i], cov = ntc_cov[i], ty = cand$type[i])
  setorder(o, k, -cov, ty)
  o[!duplicated(k), i]
}
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
  elig = Reduce(`|`, cov_arm) & ctrl_ed & cons
  gate = rep(FALSE, nrow(cand)); gate[pick(elig)] = TRUE
  gate = gate & clean[[g]]
  target = (1 - p$depletion) * ctrl$edit
  dep = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & S[[paste0("merged_", KO[i])]]$edit <= target))
  pw = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & UB[[i]] < target))
  list(gated = which(gate), true = which(gate & dep), powered = which(gate & dep & pw))
}

# ---- parity: centre at the three t equals the v9 positives --------------------------------------
v9 = fread(file.path(V9, "truth", "gt_v9_positive_sites.tsv.gz"))
key = paste(cand$chr, cand$pos, cand$strand)
par = rbindlist(lapply(GATES, function(g) rbindlist(lapply(c(t25 = 0.25, t50 = 0.50, t75 = 0.75), function(d) {
  p = modifyList(CENTRE, list(depletion = d)); tr = key[truth_of(g, p)$true]
  tn = c(`0.25` = "t25", `0.5` = "t50", `0.75` = "t75")[[as.character(d)]]
  ref = v9[gate == g & get(paste0("class_", tn)) == "positive", paste(chr, pos, read_strand)]
  data.table(gate = g, t = tn, n_lattice = length(tr), n_v9 = length(ref), n_shared = length(intersect(tr, ref)))
}))))
par[, identical := n_lattice == n_v9 & n_shared == n_v9]
print(par); fwrite(par, file.path(OUT, "parity.tsv"), sep = "\t")
stopifnot(par$identical)
msg("PARITY OK: centre reproduces the v9 positives at t25 / t50 / t75, both gates")

# ---- centre neighbours ---------------------------------------------------------------------------
nb_all = list()
for (g in GATES) {
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
    " (m = 4.6: ", mc[abs(ko_min_expected - round(-log(0.01), 4)) < 1e-9, n_true],
    "); KO cov curve: TRUE at 10: ", kc[ko_min_cov == 10, n_true])
msg("done -> ", OUT)
