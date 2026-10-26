#!/usr/bin/env Rscript
# =============================================================================
# gt_spectrum12_eval_v6.R
#
# The frozen v6 ground-truth rules applied to each of the 12 mismatch types,
# for every depletion threshold d of the step-6 sweep (0.25-0.90 by 0.05).
# Driver: human/pm/scripts/run_gt_spectrum12_v6.sh (stage eval).
#
# Per type R>Y (RNA strand), a candidate is a (position, strand) from the scan
# of the merged NTC BAM (v6 floor, per type), dbSNP common positions removed.
# Counts per unit follow parse_pileup_stranded.pl generalised to R>Y:
#   '+': genomic ref R,          cov = '.' + Y,              alt = Y
#   '-': genomic ref comp(R),    cov = ',' + lower(comp(Y)), alt = lower(comp(Y))
#   edit = alt / cov, rounded to 4 decimals as in the v6 pileups (0 if cov = 0)
# Rules (gt_robustness_grid_v6.R at the frozen centre, FROZEN_RULES.tsv):
#   covered(arm) = NTC cov >= 10 & arm cov >= 0 & arm cov x NTC edit >= 3, arms ADAR1_KO, ADARB1_KO
#   ctrl_edited  = NTC alt >= 3 & NTC edit >= 0.05
#   per teratoma (T1-T4 NTC): evaluable = cov >= 10; positive = evaluable & alt >= 3 & edit >= 0.05
#   consistent   = n_pos >= 2 & n_pos >= max(2, n_eval - 1)
#   gate         = any covered arm & ctrl_edited & consistent
#   depleted(d)  = any covered arm with KO edit <= (1 - d) x NTC edit
#   TRUE         = gate & depleted(d);  not depleted = gate & !depleted(d)
#   informative negative (2026-09-29, GROUND_TRUTH_DEFINITION.md s4) = gate & in BOTH arms
#                  covered and CP lower bound (alpha 0.05) on KO edit > 0.75 x NTC edit
# Parity: type AG at d = 0.50 must equal the v6 TRUE set (truth/gt_truth_sites.tsv.gz), and
#   AG informative negatives must equal scoring/sweep_infneg/truth_sweep_classes.tsv.gz.
#
# Outputs (gt_eval_v6/spectrum12/eval/):
#   spectrum12_by_d.tsv       type, d, n_candidates, n_gate, n_true, n_not_depleted, n_informative_negative
#   truth12_sites.tsv.gz      every gate site with its per-arm editing / KO counts, the largest d it passes,
#                             and the informative-negative flag
#   parity_AG_d050.tsv, parity_AG_informative_negative.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
V6 = file.path(B, "human/pm/gt_eval_v6")
SP = Sys.getenv("GT_SP12_DIR", file.path(V6, "spectrum12"))   # split-by-cells test: gt_eval_v8/split_cells/spectrum12
# GT_COVER_GATE=expected (default; the rules above) or depth (2026-10-01, GROUND_TRUTH_DEFINITION.md s12):
#   covered(arm) = NTC cov >= 10 & arm cov >= 10, no expected-read rule; outputs to spectrum12/eval_depth/,
#   AG parity against scoring/sweep_infneg_depth/truth_sweep_classes.tsv.gz.
#   depth<N> (2026-10-17): the depth rule with KO cov >= N; outputs to spectrum12/eval_depth<N>/, no sweep parity.
COVER_GATE = Sys.getenv("GT_COVER_GATE", "expected")
stopifnot(COVER_GATE == "expected" || grepl("^depth[0-9]*$", COVER_GATE))
MIN_COV_KO_DEPTH = if (COVER_GATE %in% c("expected", "depth")) 10 else as.numeric(sub("^depth", "", COVER_GATE))
OUT = file.path(SP, if (COVER_GATE == "expected") "eval" else paste0("eval_", COVER_GATE))
SWEEP_CLASSES = if (COVER_GATE %in% c("expected", "depth"))
  file.path(V6, "scoring", if (COVER_GATE == "expected") "sweep_infneg" else "sweep_infneg_depth",
            "truth_sweep_classes.tsv.gz") else NA_character_
dir.create(OUT, showWarnings = FALSE)
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

fr = fread(file.path(V6, "FROZEN_RULES.tsv"))
P = setNames(fr$value, fr$parameter)
MIN_COV = as.numeric(P[["min_cov"]]); MIN_COV_KO = as.numeric(P[["min_cov_ko"]])
M_EXP = as.numeric(P[["ko_min_expected"]]); K = as.numeric(P[["min_edit_reads"]])
E = as.numeric(P[["min_ctrl_edit"]]); TCOV = as.numeric(P[["min_cov_teratoma"]])
stopifnot(P[["consistency_mode"]] == "evaluable", P[["snp_filter"]] == "excluded", P[["bq"]] == "30")
D_SWEEP = round(seq(0.25, 0.90, by = 0.05), 2)
KO = c("ADAR1_KO", "ADARB1_KO"); TT = c("T1", "T2", "T3", "T4")
UNITS = c("merged_NTC_ctrl", paste0("merged_", KO), paste0(TT, "_NTC_ctrl"))
COMP = c(A = "T", C = "G", G = "C", T = "A")
BASES = c("A", "C", "G", "T")
TYPES = CJ(R = BASES, Y = BASES)[R != Y]
TYPES[, type := paste0(R, Y)]

# ---- position index (non-dbSNP union) --------------------------------------
pos = fread(file.path(SP, "positions_noSNP.bed"), header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
pos[, idx := .I]
setkey(pos, chr, pos)
N = nrow(pos); msg("positions: ", N)

# ---- raw counts per unit, aligned to pos -------------------------------------
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

# ---- candidates per type ------------------------------------------------------
cand = fread(cmd = paste("gzip -cd", shQuote(file.path(SP, "candidates_all.tsv.gz"))), sep = "\t", header = FALSE,
             select = 1:5, col.names = c("chr", "pos", "strand", "rna_ref", "rna_alt"), showProgress = FALSE)
cand = pos[cand, on = .(chr, pos), nomatch = NULL]
cand[, type := paste0(rna_ref, rna_alt)]
msg("candidate (position, strand, alt) rows after dbSNP removal: ", nrow(cand))

counts_for = function(u, idx, strand, R, Y) {
  m = U[[u]]$m
  plus = strand == "+"
  alt = ifelse(plus, m[cbind(idx, match(Y, cols))], m[cbind(idx, match(tolower(COMP[[Y]]), cols))])
  refn = ifelse(plus, m[idx, "ref_plus"], m[idx, "ref_minus"])
  cov = refn + alt
  edit = ifelse(cov > 0, as.numeric(sprintf("%.4f", alt / pmax(cov, 1))), 0)
  list(cov = cov, alt = alt, edit = edit)
}

NEG_MAX_DEP = 0.25; CI_ALPHA = 0.05
cp_lower = function(k, n, a = CI_ALPHA) { o = rep(0, length(k)); ok = n > 0 & k > 0; o[ok] = qbeta(a, k[ok], n[ok] - k[ok] + 1); o }
res = list(); gate_sites = list(); parity = NULL; parity_neg = NULL; parity_cls = list()
for (ti in seq_len(nrow(TYPES))) {
  R = TYPES$R[ti]; Y = TYPES$Y[ti]; ty = TYPES$type[ti]
  cd = cand[type == ty]
  if (!nrow(cd)) next
  gref = ifelse(cd$strand == "+", R, COMP[[R]])
  stopifnot(all(REF[cd$idx] == gref, na.rm = TRUE))
  S = lapply(setNames(UNITS, UNITS), function(u) counts_for(u, cd$idx, cd$strand, R, Y))
  ctrl = S[["merged_NTC_ctrl"]]
  cov_arm = lapply(KO, function(a) { k = S[[paste0("merged_", a)]]
    if (COVER_GATE == "expected") ctrl$cov >= MIN_COV & k$cov >= MIN_COV_KO & k$cov * ctrl$edit >= M_EXP
    else ctrl$cov >= MIN_COV & k$cov >= MIN_COV_KO_DEPTH })
  evaluable = Reduce(`|`, cov_arm)
  ctrl_ed = ctrl$alt >= K & ctrl$edit >= E
  n_eval = Reduce(`+`, lapply(TT, function(t) S[[paste0(t, "_NTC_ctrl")]]$cov >= TCOV))
  n_pos = Reduce(`+`, lapply(TT, function(t) { x = S[[paste0(t, "_NTC_ctrl")]]; x$cov >= TCOV & x$alt >= K & x$edit >= E }))
  cons = n_pos >= 2 & n_pos >= pmax(2, n_eval - 1)
  gate = evaluable & ctrl_ed & cons
  # informative negative (GROUND_TRUTH_DEFINITION.md s4): covered in both arms and the one-sided
  # Clopper-Pearson lower bound on KO editing > (1 - NEG_MAX_DEP) x NTC editing in both
  nd = lapply(seq_along(KO), function(i) { k = S[[paste0("merged_", KO[i])]]
    cov_arm[[i]] & cp_lower(k$alt, k$cov) > (1 - NEG_MAX_DEP) * ctrl$edit })
  neg = gate & nd[[1]] & nd[[2]]
  # largest d each gate site passes: max over covered arms of 1 - KO edit / NTC edit
  dmax = Reduce(pmax, lapply(seq_along(KO), function(i) {
    k = S[[paste0("merged_", KO[i])]]
    ifelse(cov_arm[[i]], 1 - k$edit / ctrl$edit, -Inf)
  }))
  for (d in D_SWEEP) {
    target = (1 - d) * ctrl$edit
    dep = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & S[[paste0("merged_", KO[i])]]$edit <= target))
    res[[length(res) + 1]] = data.table(type = ty, d = d, n_candidates = nrow(cd), n_gate = sum(gate),
                                        n_true = sum(gate & dep), n_not_depleted = sum(gate & !dep),
                                        n_informative_negative = sum(neg))
    if (ty == "AG" && !is.na(SWEEP_CLASSES) && any(abs(d - c(0.25, 0.5, 0.75)) < 1e-9)) {
      sw = fread(cmd = paste("gzip -cd", shQuote(SWEEP_CLASSES)))
      cl = paste0("class_t", round(100 * d))
      ref_set = sw[get(cl) == "positive", site]
      mine = paste0(cd$chr, ":", cd$pos)[gate & dep]
      parity_cls[[cl]] = data.table(set = paste0("positive_", cl), n_sweep = length(ref_set), n_spectrum12 = length(mine),
                                    n_shared = length(intersect(ref_set, mine)))
      if (abs(d - 0.25) < 1e-9) {
        el = sw[eligible == TRUE, site]; mg = paste0(cd$chr, ":", cd$pos)[gate]
        parity_cls[["eligible"]] = data.table(set = "eligible", n_sweep = length(el), n_spectrum12 = length(mg),
                                              n_shared = length(intersect(el, mg)))
      }
    }
    if (COVER_GATE == "expected" && ty == "AG" && abs(d - 0.5) < 1e-9) {
      v6 = fread(cmd = paste("gzip -cd", shQuote(file.path(V6, "truth", "gt_truth_sites.tsv.gz"))),
                 select = c("site", "class"))[class == "TRUE", site]
      mine = paste0(cd$chr, ":", cd$pos)[gate & dep]
      parity = data.table(n_v6 = length(v6), n_spectrum12 = length(mine), n_shared = length(intersect(v6, mine)),
                          v6_only = length(setdiff(v6, mine)), spectrum12_only = length(setdiff(mine, v6)))
      msg("PARITY AG d=0.50: ", paste(names(parity), unlist(parity), sep = "=", collapse = " "))
    }
  }
  g = which(gate)
  gate_sites[[ty]] = data.table(type = ty, chr = cd$chr[g], pos = cd$pos[g], strand = cd$strand[g],
                                ntc_cov = ctrl$cov[g], ntc_alt = ctrl$alt[g], ntc_edit = ctrl$edit[g],
                                adar1_edit = S$merged_ADAR1_KO$edit[g], adarb1_edit = S$merged_ADARB1_KO$edit[g],
                                adar1_cov = S$merged_ADAR1_KO$cov[g], adar1_alt = S$merged_ADAR1_KO$alt[g],
                                adarb1_cov = S$merged_ADARB1_KO$cov[g], adarb1_alt = S$merged_ADARB1_KO$alt[g],
                                d_max = dmax[g], informative_negative = neg[g])
  if (ty == "AG" && !is.na(SWEEP_CLASSES)) {
    sw = fread(cmd = paste("gzip -cd", shQuote(SWEEP_CLASSES)),
               select = c("site", "informative_negative"))[informative_negative == TRUE, site]
    mine = paste0(cd$chr, ":", cd$pos)[neg]
    parity_neg = data.table(n_sweep = length(sw), n_spectrum12 = length(mine), n_shared = length(intersect(sw, mine)))
    msg("PARITY AG informative negatives: ", paste(names(parity_neg), unlist(parity_neg), sep = "=", collapse = " "))
  }
  msg(ty, ": candidates ", nrow(cd), ", gate ", sum(gate), ", TRUE at d=0.5 ",
      res[[length(res) - length(D_SWEEP) + which(abs(D_SWEEP - 0.5) < 1e-9)]]$n_true)
}
out = rbindlist(res)
fwrite(out, file.path(OUT, "spectrum12_by_d.tsv"), sep = "\t")
fwrite(rbindlist(gate_sites), file.path(OUT, "truth12_sites.tsv.gz"), sep = "\t")
if (!is.null(parity)) fwrite(parity, file.path(OUT, "parity_AG_d050.tsv"), sep = "\t")
if (!is.null(parity_neg)) fwrite(parity_neg, file.path(OUT, "parity_AG_informative_negative.tsv"), sep = "\t")
pc = rbindlist(parity_cls)
fwrite(pc, file.path(OUT, "parity_AG_sweep_classes.tsv"), sep = "\t")
print(pc)
if (!is.na(SWEEP_CLASSES) && (nrow(pc) != 4 || any(pc$n_sweep != pc$n_spectrum12 | pc$n_sweep != pc$n_shared)))
  msg("WARNING: AG classes do not reproduce ", SWEEP_CLASSES, "; see parity_AG_sweep_classes.tsv")
if (COVER_GATE == "expected" && (is.null(parity) || parity$v6_only > 0 || parity$spectrum12_only > 0))
  msg("WARNING: AG at d = 0.50 does not reproduce the v6 truth set exactly; see parity_AG_d050.tsv")
msg("done")
