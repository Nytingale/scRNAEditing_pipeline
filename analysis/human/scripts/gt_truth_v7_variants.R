#!/usr/bin/env Rscript
# =============================================================================
# gt_truth_v7_variants.R   (2026-10-10)
#
# v7 truth (RNA A>G in read orientation; gt_truth_v7.R) cleaned by the artifact-prone mask and a
# truth-level nearby-mismatch filter. Site criteria, multiallelic removal, informative negatives
# and gates are v7 (gt_eval_v7/GROUND_TRUTH_DEFINITION.md); what changes is which positions stay.
#   v7            as defined
#   v7_mask       artifact-prone positions removed: RepeatMasker satellite, simple repeat or
#                 low complexity, or within 2 Mb of a centromere (acen) or heterochromatin (gvar)
#                 band (GTDEF_artifact_null_plot.R context_of)
#   v7_nomm       downstream_filter.R STEP 7 pass 1 applied to the truth: a position is removed if
#                 >= 2 mismatch types sit within w1 bp on its read strand, counting every 12-type
#                 gate site of the same gate (truth12_sites; read orientation). w1 is searched over
#                 20-200 by 20 per gate and chosen by the pipeline rule (frac_weight = 1: highest
#                 A>G share of the 12-type t50 positives that survive; first on ties).
#   v7_mask_nomm  both (a position must pass each)
#   <variant>_pgp1  (2026-10-11) the same after removing every position carrying a PGP1 germline
#                 SNV (pgp1_snps/PGP1_hu43860C_hg38_snvs.bed, PASS and VQLOW; donor of the iPSC line).
#                 PGP1 positions are dropped from the 12-type gate sites first, so they neither stay
#                 in the truth nor count as neighbours in the nearby-mismatch filter; w1 is searched
#                 again on the remaining sites.
# Removed positions are not scored (calls there are neither TP nor FP).
#
# Input : gt_eval_v6/agtc_test/truth_sites.tsv.gz, gt_eval_v6/spectrum12/{eval,eval_depth}/truth12_sites.tsv.gz,
#         ref_annots/{rmsk,cytoBand}.txt.gz, gt_eval_v7/truth/gt_v7_truth_sites.tsv.gz (parity),
#         pgp1_snps/PGP1_hu43860C_hg38_snvs.bed
# Output: gt_eval_v7/variants/truth/
#   <variant>/gt_v7_truth_sites.tsv.gz   same columns as gt_eval_v7/truth (gt_score_v7.R GT_V7_TRUTH)
#   truth_variant_counts.tsv             variant x gate x t: eligible, positives, negatives
#   truth_site_flags.tsv.gz              gate x position: context, masked, nomm_pass, w1, pgp1,
#                                        nomm_pass_pgp1, w1_pgp1
#   truth12_site_flags.tsv.gz            12-type gate sites (read strand/type): d_max, informative_negative,
#                                        context, masked, nomm_pass, pgp1, nomm_pass_pgp1
#   nomm_w1_grid.tsv                     w1 search per gate (pgp1 = FALSE / TRUE)
#   spectrum12_cleaning.tsv              12-type positives per type x gate x t x cleaning (none /
#                                        mask / nomm / mask_nomm, each also _pgp1), for the purity comparison
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
REF = file.path(B, "human/ref_annots")
OUT = Sys.getenv("GT_V7_VARIANTS_DIR", file.path(PM, "gt_eval_v7/variants/truth")); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
AGTC = Sys.getenv("GT_AGTC_DIR", file.path(PM, "gt_eval_v6/agtc_test"))
SP12 = Sys.getenv("GT_SP12_DIR", file.path(PM, "gt_eval_v6/spectrum12"))
V7_TRUTH = Sys.getenv("GT_V7_TRUTH_DIR", file.path(PM, "gt_eval_v7/truth"))
# GT_W1_FIXED (split-by-cells test): use this w1 instead of the search (frozen from the full-data truth)
W1_FIXED = as.numeric(Sys.getenv("GT_W1_FIXED", NA))
DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
GATES = c(expected = "eval", depth = "eval_depth")
# GT_GATES=<comma list> (2026-10-17): e.g. depth5,depth3 -> spectrum12/eval_depth5, eval_depth3
if (nzchar(Sys.getenv("GT_GATES"))) {
  gl = strsplit(Sys.getenv("GT_GATES"), ",", fixed = TRUE)[[1]]
  GATES = setNames(ifelse(gl == "expected", "eval", paste0("eval_", gl)), gl)
}
W_ON = seq(20, 200, 20)
CEN_PAD = 2e6
MASK = c("Satellite", "Simple repeat / low complexity", "Pericentromeric")
CTX = c(MASK, "Alu", "Other repeat", "Unique")
CHR = paste0("chr", c(1:22, "X", "Y"))
CLEAN0 = c("none", "mask", "nomm", "mask_nomm")
CLEAN = c(CLEAN0, paste0(CLEAN0, "_pgp1"))
VAR_NAME = setNames(sub("_none", "", paste0("v7_", CLEAN)), CLEAN)
PGP1_BED = file.path(PM, "pgp1_snps", "PGP1_hu43860C_hg38_snvs.bed")
pgp1 = unique(fread(PGP1_BED, header = FALSE, select = c(1, 3), col.names = c("chr", "pos")))
pgp1[, pgp1 := TRUE]
# positions kept by cleaning `cl` (rows of a table with masked, nomm_pass, nomm_pass_pgp1, pgp1)
keep_of = function(d, cl) {
  pg = grepl("_pgp1$", cl); np = if (pg) d$nomm_pass_pgp1 else d$nomm_pass
  k = switch(sub("_pgp1$", "", cl), none = rep(TRUE, nrow(d)), mask = !d$masked, nomm = np, mask_nomm = !d$masked & np)
  if (pg) k & !d$pgp1 else k
}

# ---- artifact-prone context (GTDEF_artifact_null_plot.R) --------------------------------------
rm = fread(cmd = paste("gzip -cd", shQuote(file.path(REF, "rmsk.txt.gz"))), select = c(6, 7, 8, 12, 13),
           header = FALSE, col.names = c("chr", "start", "end", "cls", "fam"))[chr %in% CHR]
rm[, ctx := fifelse(cls == "Satellite", 1L, fifelse(cls %in% c("Simple_repeat", "Low_complexity"), 2L,
              fifelse(fam == "Alu", 4L, 5L)))]
rm[, start := start + 1L]
setkey(rm, chr, start, end)
cy = fread(cmd = paste("gzip -cd", shQuote(file.path(REF, "cytoBand.txt.gz"))), header = FALSE,
           col.names = c("chr", "start", "end", "band", "stain"))
cen = cy[stain %in% c("acen", "gvar") & chr %in% CHR,
         .(start = max(0, min(start) - CEN_PAD) + 1, end = max(end) + CEN_PAD), by = chr]
setkey(cen, chr, start, end)
context_of = function(chr, pos) {
  q = unique(data.table(chr = chr, start = pos, end = pos)); setkey(q, chr, start, end)
  r = foverlaps(q, rm, nomatch = NULL)[, .(rep = min(ctx)), by = .(chr, pos = i.start)]
  c3 = foverlaps(q, cen, nomatch = NULL)[, .(chr, pos = i.start, peri = TRUE)]
  q = merge(q[, .(chr, pos = start)], r, by = c("chr", "pos"), all.x = TRUE)
  q = merge(q, unique(c3), by = c("chr", "pos"), all.x = TRUE)
  q[, k := fifelse(!is.na(rep) & rep <= 2L, rep, fifelse(!is.na(peri), 3L, fifelse(!is.na(rep), rep, 6L)))]
  q[, .(chr, pos, context = factor(CTX[k], CTX))]
}

# ---- nearby-mismatch pass 1 against the 12-type gate sites ------------------------------------
# number of mismatch types among `nb` within +-w on the query's strand (filReliable rule)
n_types_near = function(q, nb, w) {
  q = copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h = nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL,
         .(row = i.row, type = x.type)]
  h = rbind(h, q[, .(row, type)])
  h[, .(n_types = uniqueN(type)), by = row][order(row), n_types]
}

ag = fread(file.path(AGTC, "truth_sites.tsv.gz"))
ag = ag[variant == "AG" & eligible == TRUE & multiallelic == FALSE]
ag[, c("chr", "pos") := tstrsplit(site, ":", fixed = TRUE)][, pos := as.integer(pos)]
stopifnot(all(ag$type == "AG"), all(ag$strand %in% c("+", "-")))

flags = list(); grid_all = list(); spec = list(); flags12 = list()
for (g in names(GATES)) {
  t12 = fread(file.path(SP12, GATES[[g]], "truth12_sites.tsv.gz"))
  t12 = merge(t12, context_of(t12$chr, t12$pos), by = c("chr", "pos"))
  t12[, masked := context %in% MASK]
  t12 = merge(t12, pgp1, by = c("chr", "pos"), all.x = TRUE)
  t12[is.na(pgp1), pgp1 := FALSE]

  # nearby-mismatch pass on all gate sites, then on the gate sites without PGP1 positions
  w1s = c()
  for (pg in c(FALSE, TRUE)) {
    u = if (pg) t12[pgp1 == FALSE] else copy(t12)
    nb = u[, .(chr, pos, strand, type)]
    # w1 search on the 12-type t50 positives (pipeline rule, frac_weight = 1)
    p50 = u[d_max >= DS[["t50"]], .(chr, pos, strand, type)]
    gr = rbindlist(lapply(W_ON, function(w) {
      keep = n_types_near(p50, nb, w) < 2
      data.table(gate = g, pgp1 = pg, w1 = w, n_kept = sum(keep), n_ag = sum(keep & p50$type == "AG"))
    }))
    gr[, frac_ag := n_ag / pmax(n_kept, 1)]
    sc = if (uniqueN(gr$frac_ag) == 1) rep(0, nrow(gr)) else as.numeric(scale(gr$frac_ag))
    w1_best = if (is.na(W1_FIXED)) gr$w1[which.max(sc)] else W1_FIXED
    gr[, chosen := w1 == w1_best]
    grid_all[[paste(g, pg)]] = gr
    message(g, if (pg) " (PGP1 removed)" else "", ": w1 = ", w1_best, "; PGP1 gate sites ", t12[, sum(pgp1)])
    ps = n_types_near(u[, .(chr, pos, strand, type)], nb, w1_best) < 2
    if (pg) {
      t12[, nomm_pass_pgp1 := FALSE]
      t12[pgp1 == FALSE, nomm_pass_pgp1 := ps]
    } else t12[, nomm_pass := ps]
    w1s[[if (pg) "pgp1" else "all"]] = w1_best
  }

  # every eligible A>G truth row of this gate must be a 12-type gate site
  tg = merge(ag[gate == g], t12[, .(chr, pos, strand, type, context, masked, nomm_pass, pgp1, nomm_pass_pgp1)],
             by = c("chr", "pos", "strand", "type"), all.x = TRUE)
  stopifnot(!anyNA(tg$masked))
  tg[, `:=`(w1 = w1s[["all"]], w1_pgp1 = w1s[["pgp1"]])]
  flags[[g]] = tg
  flags12[[g]] = t12[, .(gate = g, chr, pos, strand, type, d_max, informative_negative, context, masked,
                         nomm_pass, pgp1, nomm_pass_pgp1)]

  # 12-type positives under each cleaning
  for (tn in names(DS)) {
    s = t12[d_max >= DS[[tn]]]
    for (cl in CLEAN) {
      w = if (grepl("_pgp1$", cl)) w1s[["pgp1"]] else w1s[["all"]]
      spec[[length(spec) + 1]] = s[keep_of(s, cl), .N, by = type][, `:=`(gate = g, t = tn, cleaning = cl, w1 = w)]
    }
  }
}
fl = rbindlist(flags)
fwrite(rbindlist(grid_all), file.path(OUT, "nomm_w1_grid.tsv"), sep = "\t")
fwrite(fl[, .(gate, site, chr, pos, read_type = type, read_strand = strand, context, masked,
              nomm_pass, w1, pgp1, nomm_pass_pgp1, w1_pgp1, pos_t25, pos_t50, pos_t75, informative_negative)],
       file.path(OUT, "truth_site_flags.tsv.gz"), sep = "\t")
fwrite(rbindlist(flags12), file.path(OUT, "truth12_site_flags.tsv.gz"), sep = "\t")
sp = rbindlist(spec)
sp = merge(CJ(gate = names(GATES), t = names(DS), cleaning = CLEAN, type = sort(unique(sp$type))), sp,
           by = c("gate", "t", "cleaning", "type"), all.x = TRUE)
sp[is.na(N), N := 0L][, w1 := max(w1, na.rm = TRUE), by = .(gate, cleaning)]
sp[, pct := 100 * N / sum(N), by = .(gate, t, cleaning)]
fwrite(sp, file.path(OUT, "spectrum12_cleaning.tsv"), sep = "\t")

# ---- write the eight truth tables ----------------------------------------------------------------
fl[, genomic_ref := fifelse(strand == "+", "A", "T")]
fl[, a_to_i_strand := strand]
for (tn in names(DS)) {
  p = fl[[paste0("pos_", tn)]]
  set(fl, j = paste0("class_", tn), value = fifelse(p, "positive", fifelse(fl$informative_negative, "informative_negative", "unresolved")))
}
COLS = c("gate", "site", "chr", "pos", "genomic_ref", "a_to_i_strand", "read_type", "read_strand", paste0("class_", names(DS)))
setnames(fl, c("type", "strand"), c("read_type", "read_strand"))

v7 = fread(file.path(V7_TRUTH, "gt_v7_truth_sites.tsv.gz"))
stopifnot(all(v7$read_type == "AG"))
cnt = list()
for (cl in CLEAN) {
  nm = VAR_NAME[[cl]]
  x = fl[keep_of(fl, cl), ..COLS]; setorder(x, gate, chr, pos)
  if (cl == "none") {
    v = copy(v7); setorder(v, gate, chr, pos)
    stopifnot(isTRUE(all.equal(as.data.frame(x), as.data.frame(v[, ..COLS]), check.attributes = FALSE)))
    message("v7 parity with gt_eval_v7/truth: identical")
  }
  dir.create(file.path(OUT, nm), showWarnings = FALSE)
  fwrite(x, file.path(OUT, nm, "gt_v7_truth_sites.tsv.gz"), sep = "\t")
  cnt[[nm]] = rbindlist(lapply(names(DS), function(tn) {
    cc = paste0("class_", tn)
    x[, .(variant = nm, cleaning = cl, t = tn, eligible = .N,
          positive = sum(get(cc) == "positive"),
          informative_negative = sum(get(cc) == "informative_negative")), by = gate]
  }))
}
cnt = rbindlist(cnt)
fwrite(cnt, file.path(OUT, "truth_variant_counts.tsv"), sep = "\t")
print(dcast(cnt, gate + t ~ variant, value.var = "positive"))
print(dcast(cnt[t == "t50"], gate ~ variant, value.var = "informative_negative"))
print(dcast(sp[t == "t50", .(ag_pct = sum(pct[type == "AG"]), n_ag = sum(N[type == "AG"]), n = sum(N)),
               by = .(gate, cleaning)], gate ~ cleaning, value.var = c("ag_pct", "n_ag", "n")))
cat("wrote", OUT, "\n")
