#!/usr/bin/env Rscript
# =============================================================================
# gt_robustness_grid_v6.R
#
# GT v6 robustness analysis: is the ground-truth set robust to small
# perturbations of every threshold, while keeping a substantial number of sites?
# Evaluates the v6 truth rules (gt_eval_classify_v6.R, reproduced exactly) for
# every combination on a 7-axis lattice around a fixed centre, then compares each
# combination with its 1- and 2-step neighbours.
#
#   axis              values (centre in [])
#   min_cov           8 9 [10] 11 12      merged NTC coverage
#   min_cov_teratoma  8 9 [10] 11 12      per-teratoma (per-sample) NTC coverage, reproducibility
#   ko_min_expected   1 2 [3] 4 5         informative KO: KO cov x NTC editing (either arm)
#   min_edit_reads    2 [3] 4             NTC edited reads
#   min_ctrl_edit     .03 .04 [.05] .06 .07   NTC editing level
#   depletion         .3 .4 [.5] .6 .7    KO <= (1-d) x control, either arm
#   rule              strict2 [evaluable] strict3
#   = 28,125 combinations (5 x 5 x 5 x 3 x 5 x 5 x 3)
#   fixed: min_cov_ko 0 (no raw KO floor)
#
# 2026-09-27 (informative-KO revision): the raw KO-coverage axis (min_cov_ko
# 3-7) is replaced by ko_min_expected, the edited reads a KO arm would carry
# under no depletion. It scales with each site's editing level, so it does not
# encode this dataset's (shallow) ADAR1_KO depth.
# The earlier 7-axis KO-coverage outputs are in robustness/_superseded_kocov/.
# 2026-09-27: per-teratoma NTC coverage swept like NTC coverage (8-12, step 1;
# user: a depth-invariant version is not appropriate here). The 6-axis run with
# it fixed at 10 is in robustness/_superseded_6axis/.
#
# Rules (identical to the engine; dbSNP common removed, BQ30 pileups):
#   covered(arm) = NTC cov >= min_cov & arm cov >= min_cov_ko & arm cov x NTC editing >= m
#   evaluable    = any arm covered
#   ctrl_edited  = NTC edited reads >= k & NTC editing >= e
#   per teratoma: evaluable = cov >= tcov; positive = cov >= tcov & alt >= k & edit >= e
#   consistent   = strict m: n_positive >= m | evaluable: n_pos >= 2 & n_pos >= max(2, n_eval - 1)
#   depleted     = any covered arm with KO editing <= (1 - d) x NTC editing
#   TRUE         = evaluable & ctrl_edited & consistent & depleted
#
# Modes (driver: human/pm/scripts/run_gt_truth_rebuild_v6.sh robust*):
#   task   <i>   evaluate block i (one min_cov x ko_min_expected pair, 1,125 combos)
#   gather       neighbour Jaccard / fraction lost / site counts for every combo
#   parity       compare TRUE sets with real engine runs in robustness/parity/
#   mcurve       TRUE sites vs a fine informative-KO (m) grid at the centre -> m_curve.tsv
#
# Inputs : human/pm/gt_eval_v6/pileups/bq30/pileup_{merged_*,T*_NTC_ctrl}.tsv,
#          candidates/candidates_noSNP.sites.txt, candidates/rediportal_sites.txt
# Outputs: human/pm/gt_eval_v6/robustness/
#          blocks/block_<i>.rds (TRUE site indices per combo), blocks/block_<i>.tsv,
#          combos.tsv.gz, neighbours.tsv.gz, parity.tsv
# =============================================================================

suppressPackageStartupMessages(library(data.table))
args = commandArgs(trailingOnly = TRUE)
MODE = args[1]
B = Sys.getenv("SCRNA_DIR", ".")
V6 = file.path(B, "human/pm/gt_eval_v6")
OUT = file.path(V6, "robustness")
dir.create(file.path(OUT, "blocks"), recursive = TRUE, showWarnings = FALSE)
PILE = file.path(V6, "pileups", "bq30")
CTRL = "NTC_ctrl"; KO = c("ADAR1_KO", "ADARB1_KO"); TT = c("T1", "T2", "T3", "T4")
ALPHA = 0.05
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

AXES = list(
  min_cov          = c(8, 9, 10, 11, 12),
  min_cov_teratoma = c(8, 9, 10, 11, 12),
  ko_min_expected  = c(1, 2, 3, 4, 5),
  min_edit_reads   = c(2, 3, 4),
  min_ctrl_edit    = c(0.03, 0.04, 0.05, 0.06, 0.07),
  depletion        = c(0.3, 0.4, 0.5, 0.6, 0.7),
  rule             = c("strict2", "evaluable", "strict3"))
FIXED = list(min_cov_ko = 0)
CENTRE = list(min_cov = 10, min_cov_teratoma = 10, ko_min_expected = 3, min_edit_reads = 3, min_ctrl_edit = 0.05,
              depletion = 0.5, rule = "evaluable")
GRID = do.call(CJ, c(AXES, sorted = FALSE))
GRID[, combo := sprintf("c%05d", .I)]
GRID[, names(FIXED) := FIXED]
BLOCKS = unique(GRID[, .(min_cov, ko_min_expected)])[, block := .I]
GRID = merge(GRID, BLOCKS, by = c("min_cov", "ko_min_expected"), sort = FALSE)
setorder(GRID, combo)

# ---- grid + parity combos (for the driver) -------------------------------------
PARITY = rbindlist(list(
  as.data.table(CENTRE),
  data.table(min_cov = 12, min_cov_teratoma = 10, ko_min_expected = 1, min_edit_reads = 3, min_ctrl_edit = 0.05, depletion = 0.5, rule = "evaluable"),
  data.table(min_cov = 10, min_cov_teratoma = 8, ko_min_expected = 3, min_edit_reads = 2, min_ctrl_edit = 0.03, depletion = 0.5, rule = "evaluable"),
  data.table(min_cov = 9, min_cov_teratoma = 12, ko_min_expected = 2, min_edit_reads = 4, min_ctrl_edit = 0.07, depletion = 0.3, rule = "evaluable"),
  data.table(min_cov = 11, min_cov_teratoma = 9, ko_min_expected = 4, min_edit_reads = 3, min_ctrl_edit = 0.04, depletion = 0.6, rule = "strict2"),
  data.table(min_cov = 10, min_cov_teratoma = 11, ko_min_expected = 5, min_edit_reads = 4, min_ctrl_edit = 0.06, depletion = 0.7, rule = "strict3"),
  data.table(min_cov = 8, min_cov_teratoma = 10, ko_min_expected = 1, min_edit_reads = 2, min_ctrl_edit = 0.05, depletion = 0.4, rule = "strict2")))
if (MODE == "grid") {
  fwrite(GRID, file.path(OUT, "grid.tsv"), sep = "\t")
  pc = merge(PARITY, GRID, by = names(AXES))
  stopifnot(nrow(pc) == nrow(PARITY))
  pc[, `:=`(consistency_mode = fifelse(rule == "evaluable", "evaluable", "strict"),
            consistency_min = fifelse(rule == "evaluable", 3L, as.integer(sub("strict", "", rule))))]
  fwrite(pc[, .(combo, block, min_cov, min_cov_ko, ko_min_expected, min_edit_reads, min_ctrl_edit, depletion,
                min_cov_teratoma, consistency_mode, consistency_min)], file.path(OUT, "parity_combos.tsv"), sep = "\t")
  msg("grid: ", nrow(GRID), " combos in ", nrow(BLOCKS), " blocks; parity combos: ", nrow(pc))
}

# ---- site-level inputs, aligned to the non-dbSNP candidate universe ----------
load_sites = function() {
  sites = fread(file.path(V6, "candidates", "candidates_noSNP.sites.txt"), header = FALSE)$V1
  rd = function(unit) {
    p = fread(file.path(PILE, paste0("pileup_", unit, ".tsv")),
              select = c("chr", "pos", "coverage", "variant_cov", "editing_level"), showProgress = FALSE)
    p[, site := paste0(chr, ":", pos)]
    m = match(sites, p$site)
    list(cov = fifelse(is.na(m), 0, as.numeric(p$coverage[m])),
         alt = fifelse(is.na(m), 0, as.numeric(p$variant_cov[m])),
         edit = fifelse(is.na(m), 0, suppressWarnings(as.numeric(p$editing_level[m]))))
  }
  S = list(site = sites, ctrl = rd(paste0("merged_", CTRL)),
           ko = lapply(setNames(KO, KO), function(a) rd(paste0("merged_", a))),
           tt = lapply(setNames(TT, TT), function(t) rd(paste0(t, "_", CTRL))))
  for (a in KO) S$ko[[a]]$edit[is.na(S$ko[[a]]$edit)] = 0
  S$ctrl$edit[is.na(S$ctrl$edit)] = 0
  redi = fread(file.path(V6, "candidates", "rediportal_sites.txt"), header = FALSE)$V1
  S$redi = sites %in% redi
  S
}
cp0 = function(n) { out = rep(1, length(n)); ok = n > 0; out[ok] = qbeta(1 - ALPHA, 1, n[ok]); out }

# ---- evaluate every combo of one block ----------------------------------------
eval_block = function(S, bi) {
  g = GRID[block == bi]
  mc = g$min_cov[1]; mk = g$min_cov_ko[1]; me = g$ko_min_expected[1]
  cov_arm = lapply(KO, function(a) S$ctrl$cov >= mc & S$ko[[a]]$cov >= mk & S$ko[[a]]$cov * S$ctrl$edit >= me)
  evaluable = Reduce(`|`, cov_arm)
  ub = lapply(KO, function(a) cp0(S$ko[[a]]$cov))
  sets = vector("list", nrow(g)); names(sets) = g$combo
  summ = vector("list", nrow(g))
  for (k in AXES$min_edit_reads) for (e in AXES$min_ctrl_edit) {
    ctrl_ed = S$ctrl$alt >= k & S$ctrl$edit >= e
    base_ke = evaluable & ctrl_ed
    for (tc in AXES$min_cov_teratoma) {
      n_eval = Reduce(`+`, lapply(TT, function(t) S$tt[[t]]$cov >= tc))
      n_pos = Reduce(`+`, lapply(TT, function(t) S$tt[[t]]$cov >= tc & S$tt[[t]]$alt >= k & S$tt[[t]]$edit >= e))
      n_pos[is.na(n_pos)] = 0
      cons = list(strict2 = n_pos >= 2, strict3 = n_pos >= 3,
                  evaluable = n_pos >= 2 & n_pos >= pmax(2, n_eval - 1))
      for (d in AXES$depletion) {
        target = (1 - d) * S$ctrl$edit
        dep = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & S$ko[[KO[i]]]$edit <= target))
        pw = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & ub[[i]] < target))
        for (ru in AXES$rule) {
          j = which(g$min_edit_reads == k & abs(g$min_ctrl_edit - e) < 1e-9 & g$min_cov_teratoma == tc &
                      abs(g$depletion - d) < 1e-9 & g$rule == ru)
          gate = base_ke & cons[[ru]]
          tr = gate & dep
          sets[[j]] = which(tr)
          summ[[j]] = data.table(combo = g$combo[j], n_true = sum(tr), n_true_powered = sum(tr & pw),
                                 n_fp_known = sum(gate & !dep & S$redi), n_fp_novel = sum(gate & !dep & !S$redi),
                                 pct_true_rediportal = if (sum(tr)) 100 * mean(S$redi[tr]) else NA_real_)
        }
      }
    }
  }
  list(sets = sets, summary = merge(g, rbindlist(summ), by = "combo"))
}

if (MODE == "task") {
  bi = as.integer(args[2]); stopifnot(bi %in% BLOCKS$block)
  msg("block ", bi, ": min_cov=", BLOCKS$min_cov[bi], " ko_min_expected=", BLOCKS$ko_min_expected[bi])
  S = load_sites(); msg("sites: ", length(S$site))
  r = eval_block(S, bi)
  saveRDS(r$sets, file.path(OUT, "blocks", sprintf("block_%02d.rds", bi)))
  fwrite(r$summary, file.path(OUT, "blocks", sprintf("block_%02d.tsv", bi)), sep = "\t")
  msg("done: ", nrow(r$summary), " combos; centre-block TRUE range ", min(r$summary$n_true), "-", max(r$summary$n_true))
}

# ---- neighbours: 1- and 2-step moves along each axis ----------------------------
if (MODE == "gather") {
  files = file.path(OUT, "blocks", sprintf("block_%02d.rds", BLOCKS$block))
  miss = files[!file.exists(files)]; if (length(miss)) stop("missing blocks: ", paste(basename(miss), collapse = ", "))
  sets = do.call(c, lapply(files, readRDS))
  summ = rbindlist(lapply(sub("\\.rds$", ".tsv", files), fread))
  setkey(summ, combo)
  N = length(fread(file.path(V6, "candidates", "candidates_noSNP.sites.txt"), header = FALSE)$V1)
  pos = GRID[, lapply(names(AXES), function(a) match(get(a), AXES[[a]]))]; setnames(pos, names(AXES))
  key_of = function(p) do.call(paste, c(as.list(p), sep = "_"))
  idx = setNames(GRID$combo, key_of(pos))
  msg("gather: ", length(sets), " combos")
  nb = vector("list", nrow(GRID))
  for (i in seq_len(nrow(GRID))) {
    A = sets[[GRID$combo[i]]]; bit = logical(N); bit[A] = TRUE; nA = length(A)
    rows = list()
    for (a in names(AXES)) for (s in c(-2, -1, 1, 2)) {
      p = unlist(pos[i]); p[[a]] = p[[a]] + s
      if (p[[a]] < 1 || p[[a]] > length(AXES[[a]])) next
      cb = idx[[key_of(p)]]; Bs = sets[[cb]]; inter = sum(bit[Bs])
      rows[[length(rows) + 1]] = list(combo = GRID$combo[i], axis = a, step = s, neighbour = cb,
                                      neighbour_value = as.character(AXES[[a]][p[[a]]]),
                                      n_true = nA, n_neighbour = length(Bs),
                                      jaccard = if (nA + length(Bs) - inter > 0) inter / (nA + length(Bs) - inter) else NA_real_,
                                      frac_lost = if (nA > 0) 1 - inter / nA else NA_real_)
    }
    nb[[i]] = rbindlist(rows)
    if (i %% 2000 == 0) msg("  ", i)
  }
  nb = rbindlist(nb)
  fwrite(nb, file.path(OUT, "neighbours.tsv.gz"), sep = "\t")
  agg = nb[, .(min_jaccard_1 = min(jaccard[abs(step) == 1]), min_jaccard_2 = min(jaccard),
               max_frac_lost_1 = max(frac_lost[abs(step) == 1]), max_frac_lost_2 = max(frac_lost),
               min_neighbour_sites_1 = min(n_neighbour[abs(step) == 1]),
               n_neighbours_1 = sum(abs(step) == 1), n_neighbours_2 = .N), by = combo]
  out = merge(summ, agg, by = "combo")
  out[, `:=`(
             is_centre = min_cov == CENTRE$min_cov & min_cov_teratoma == CENTRE$min_cov_teratoma &
               ko_min_expected == CENTRE$ko_min_expected &
               min_edit_reads == CENTRE$min_edit_reads & abs(min_ctrl_edit - CENTRE$min_ctrl_edit) < 1e-9 &
               abs(depletion - CENTRE$depletion) < 1e-9 & rule == CENTRE$rule,
             # every axis has both 1-step neighbours / every possible <= 2-step neighbour
             interior_1 = n_neighbours_1 == 2 * length(AXES),
             interior_2 = n_neighbours_2 == sum(pmin(4, 2 * (lengths(AXES) - 1))))]
  fwrite(out, file.path(OUT, "combos.tsv.gz"), sep = "\t")
  c0 = out[is_centre == TRUE]
  msg("centre: n_true=", c0$n_true, " min J (1-step)=", round(c0$min_jaccard_1, 3), " min J (2-step)=",
      round(c0$min_jaccard_2, 3), " max lost (1-step)=", round(c0$max_frac_lost_1, 3))
  msg("max over combos of min 1-step J: ", round(max(out$min_jaccard_1, na.rm = TRUE), 3))
}

# ---- m curve: TRUE sites vs a fine informative-KO grid at the centre -------------
# Every other threshold at CENTRE; m on a 0.1 grid plus exactly -log(alpha) for
# alpha = 0.05 / 0.01 (m = 3.00 / 4.61), where exp(-m) = alpha.
if (MODE == "mcurve") {
  S = load_sites(); msg("sites: ", length(S$site))
  k = CENTRE$min_edit_reads; e = CENTRE$min_ctrl_edit; d = CENTRE$depletion; tc = CENTRE$min_cov_teratoma
  n_eval = Reduce(`+`, lapply(TT, function(t) S$tt[[t]]$cov >= tc))
  n_pos = Reduce(`+`, lapply(TT, function(t) S$tt[[t]]$cov >= tc & S$tt[[t]]$alt >= k & S$tt[[t]]$edit >= e))
  gate = S$ctrl$cov >= CENTRE$min_cov & S$ctrl$alt >= k & S$ctrl$edit >= e &
    n_pos >= 2 & n_pos >= pmax(2, n_eval - 1)
  target = (1 - d) * S$ctrl$edit
  ub = lapply(KO, function(a) cp0(S$ko[[a]]$cov))
  ms = sort(unique(round(c(seq(0, 6, 0.1), -log(0.05), -log(0.01)), 4)))
  res = rbindlist(lapply(ms, function(m) {
    cov_arm = lapply(KO, function(a) S$ctrl$cov >= CENTRE$min_cov & S$ko[[a]]$cov * S$ctrl$edit >= m & S$ko[[a]]$cov > 0)
    ev = Reduce(`|`, cov_arm)
    dep = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & S$ko[[KO[i]]]$edit <= target))
    pw = Reduce(`|`, lapply(seq_along(KO), function(i) cov_arm[[i]] & ub[[i]] < target))
    data.table(ko_min_expected = m, n_gated = sum(gate & ev), n_true = sum(gate & ev & dep),
               n_true_powered = sum(gate & ev & dep & pw))
  }))
  fwrite(res, file.path(OUT, "m_curve.tsv"), sep = "\t")
  msg("m curve: ", nrow(res), " values; TRUE at m=3: ", res[abs(ko_min_expected - 3) < 1e-9, n_true],
      ", m=4.6052: ", res[abs(ko_min_expected - 4.6052) < 1e-4, n_true])
}

# ---- parity with real engine runs ------------------------------------------------
# robustness/parity/<combo>/gt_truth_sites.tsv.gz written by the driver's
# `robustparity` stage (engine with that combo's parameters).
if (MODE == "parity") {
  sites = fread(file.path(V6, "candidates", "candidates_noSNP.sites.txt"), header = FALSE)$V1
  dirs = list.dirs(file.path(OUT, "parity"), recursive = FALSE)
  res = rbindlist(lapply(dirs, function(dd) {
    cb = basename(dd)
    eng = fread(file.path(dd, "gt_truth_sites.tsv.gz"), select = c("site", "class"))[class == "TRUE", site]
    bi = GRID[combo == cb, block]
    ev = sites[readRDS(file.path(OUT, "blocks", sprintf("block_%02d.rds", bi)))[[cb]]]
    data.table(combo = cb, n_engine = length(eng), n_evaluator = length(ev),
               n_shared = length(intersect(eng, ev)), identical = setequal(eng, ev))
  }))
  print(res)
  fwrite(res, file.path(OUT, "parity.tsv"), sep = "\t")
  if (!all(res$identical)) stop("evaluator != engine for: ", paste(res[identical == FALSE, combo], collapse = ", "))
  msg("PARITY OK: ", nrow(res), " combos identical to the engine")
}
