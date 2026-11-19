#!/usr/bin/env Rscript
# =============================================================================
# gt_v9_sj_indel_calls.R   (2026-11-16; v9 port of gt_v8_sj_indel_calls.R; Reviewer 2, major 3; Response Figure 17 text)
# Do filters that remove calls near splice junctions or indels raise F1? Merged NTC arm tier, the 76 ranked
# configurations, main call state (step 5), v9 truth and scoring (gt_v9_build.sh: GT_FP_RULE=ag, GT_UNIT_RULE=gate_cov).
#   d_sj     distance to the nearest splice junction boundary (first or last intron base of a STAR junction with >= 5
#            uniquely mapped reads summed over the 8 teratoma runs, or annotated), as gt_v8_sj_indel_calls.R
#   d_indel  distance to the nearest indel with >= MIN_INDEL reads in the seed-8787 merged NTC BAM
#            (sj_indel_v9.sh -> gt_eval_v9/indel_sj/indels); deletion = its first and last deleted base,
#            insertion = the two reference bases around it
# Filters: remove calls with d_sj <= d (d = 1-6, 8, 10) or d_indel <= d (d = 1, 2, 3, 5, 10); per configuration the
# A-to-G calls removed and the arm F1 for the six KO truth sets are recomputed (TP: A-to-G call at a positive
# covered in the arm; FP: A-to-G call at an informative negative). Unfiltered F1 must equal
# gt_eval_v9/scoring_pwcov76/mismatch/arm/gt_v7_unit_scores.tsv.
# Calls: <root>/NTC_ctrl/04_filtered/<caller>/<bam_type>/<strand_mode>/snp_<x>/*_step5_nomm.tsv.gz for the roots
#   results_gt_v8 (58), gt_eval_v8/cellsnp_hybrid/results/centre/results_gt_v8 (2), results_gt_v8_umimd (16)
# Outputs (human/pm/gt_eval_v9/indel_sj/):
#   filter_compare.tsv   configuration x filter x d (d = 0: no filter): calls, A-to-G calls, F1 / precision / recall
#   filter_summary.tsv   filter x d: A-to-G calls removed (min / median / max over configurations), F1 change over
#                        configurations x truth sets (median, max, number raised)
# Usage : srun ... Rscript human/scripts/gt_v9_sj_indel_calls.R   (after sj_indel_v9.sh submit)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("SCRNA_DIR", ".")
PM <- file.path(ROOT, "human/pm")
V9 <- file.path(PM, "gt_eval_v9")
SC <- file.path(V9, "scoring_pwcov76", "mismatch")
OUT <- file.path(V9, "indel_sj")
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
MIN_INDEL <- as.integer(Sys.getenv("MIN_INDEL", "3"))
MIN_COV <- 10; MIN_EXP <- 3
SJ_D <- c(1:6, 8L, 10L)
INDEL_D <- c(1L, 2L, 3L, 5L, 10L)
ROOTS <- file.path(c(file.path(PM, "results_gt_v8"), file.path(PM, "gt_eval_v8/cellsnp_hybrid/results/centre/results_gt_v8"),
                     file.path(PM, "results_gt_v8_umimd")), "NTC_ctrl", "04_filtered")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

ranked <- unique(fread(file.path(SC, "config_rank.tsv"))[, c(COMBO, "config"), with = FALSE])
stopifnot(nrow(ranked) == 76L)
meta <- rbindlist(lapply(ROOTS, function(r) {
  f <- list.files(r, recursive = TRUE, full.names = TRUE, pattern = "_step5_nomm\\.tsv\\.gz$")
  m <- data.table(f = f)
  m[, c("caller", "bam_type", "strand_mode", "snp_dir") := tstrsplit(sub(paste0(r, "/"), "", f, fixed = TRUE), "/",
                                                                      fixed = TRUE, keep = 1:4)]
  m
}))
meta[, snp_subtract := sub("snp_", "", snp_dir)]
meta[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
meta <- meta[config %in% ranked$config]
stopifnot(nrow(meta) == 76L, uniqueN(meta$config) == 76L)
calls <- unique(rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  x <- fread(meta$f[i], select = c("chr", "pos", "varType", "strand"), showProgress = FALSE)
  x[, config := meta$config[i]]
})))
calls[, ag := varType == "AG"]
msg("calls:", nrow(calls), "; A-to-G:", sum(calls$ag))

# ---- distances ------------------------------------------------------------------
nearest <- function(pts, q) {
  setkey(pts, chr, b)
  qq <- q[, .(chr, b = pos, pos)]
  abs(pts[qq, on = .(chr, b), roll = "nearest", x.b] - qq$pos)
}
sjf <- list.files(file.path(PM, "01_starsolobams"), "_SJ\\.out\\.tab$", full.names = TRUE)
stopifnot(length(sjf) == 8)
sj <- rbindlist(lapply(sjf, fread, col.names = c("chr", "start", "end", "strand", "motif", "annot", "uniq", "multi",
                                                 "overhang")))
sj <- sj[, .(uniq = sum(uniq), annot = max(annot)), by = .(chr, start, end)][uniq >= 5 | annot == 1]
bd <- unique(rbind(sj[, .(chr, b = start)], sj[, .(chr, b = end)]))
upos <- unique(calls[, .(chr, pos)])
upos[, d_sj := nearest(bd, upos)]
inf <- list.files(file.path(OUT, "indels"), "\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(inf) == 23)
ind <- rbindlist(lapply(inf, fread, col.names = c("chr", "pos", "type", "len", "n_reads")))
ind <- ind[, .(n_reads = sum(n_reads)), by = .(chr, pos, type, len)][n_reads >= MIN_INDEL]
ib <- unique(rbind(ind[type == "D", .(chr, b = pos)], ind[type == "D", .(chr, b = pos + len - 1L)],
                   ind[type == "I", .(chr, b = pos - 1L)], ind[type == "I", .(chr, b = pos)]))
upos[, d_indel := nearest(ib, upos)]
upos[, site := paste0(chr, ":", pos)]
calls <- upos[calls, on = .(chr, pos)]

# ---- arm F1 (gt_score_v7.R, GT_FP_RULE=ag, GT_UNIT_RULE=gate_cov) -----------------------
truth <- fread(file.path(V9, "truth", "gt_v9_truth_sites.tsv.gz"), select = c("gate", "site", paste0("class_t", c(25, 50, 75))))
stopifnot(!anyDuplicated(truth[, .(gate, site)]))
tl <- melt(truth, id.vars = c("gate", "site"), variable.name = "t", value.name = "cls")[cls %chin% c("positive", "informative_negative")]
tl[, t := sub("class_", "", t)]
uc <- fread(file.path(V9, "unit_cov", "ko_truth_unit_cov.tsv.gz"))[unit == "NTC_ctrl"]
covered <- rbind(uc[coverage >= MIN_COV, .(site, gate = "depth")],
                 uc[coverage >= MIN_COV & coverage > 0 & coverage * (variant_cov / coverage) >= MIN_EXP, .(site, gate = "expected")])
tl[, covered := FALSE][covered, covered := TRUE, on = .(site, gate)]
ref <- fread(file.path(SC, "arm", "gt_v7_unit_scores.tsv"))[unit == "NTC_ctrl"]
npu <- unique(ref[, .(gate, t, n_pos_unit)])
stopifnot(!anyDuplicated(npu[, .(gate, t)]))
ct <- unique(calls[site %chin% truth$site & ag == TRUE, .(config, site, d_sj, d_indel)])
f1_of <- function(k) {
  x <- merge(unique(k[, .(config, site)]), tl, by = "site", allow.cartesian = TRUE)
  s <- x[, .(n_tp = sum(cls == "positive" & covered), n_fp = sum(cls == "informative_negative")), by = .(config, gate, t)]
  s <- merge(CJ(config = ranked$config, gate = c("depth", "expected"), t = c("t25", "t50", "t75")), s,
             by = c("config", "gate", "t"), all.x = TRUE)
  for (v in c("n_tp", "n_fp")) set(s, which(is.na(s[[v]])), v, 0L)
  s <- merge(s, npu, by = c("gate", "t"))
  s[, `:=`(precision = fifelse(n_tp + n_fp > 0, n_tp / (n_tp + n_fp), NA_real_), recall = n_tp / n_pos_unit)]
  s[, f1 := fifelse(is.na(precision) | precision + recall == 0, 0, 2 * precision * recall / (precision + recall))]
  s
}
base <- f1_of(ct)
chk <- merge(base, ref[, .(config, gate, t, f1_ref = f1, tp_ref = n_tp_cov, fp_ref = n_fp)], by = c("config", "gate", "t"))
stopifnot(nrow(chk) == 76L * 6L)
msg("arm F1 vs v9 scoring: max |diff| =", signif(chk[, max(abs(f1 - f1_ref))], 3),
    "; TP diff", chk[, sum(n_tp != tp_ref)], "; FP diff", chk[, sum(n_fp != fp_ref)])
stopifnot(chk[, max(abs(f1 - f1_ref))] < 1e-9)

# ---- filters --------------------------------------------------------------------------
grid <- rbind(data.table(filter = "none", d = 0L), data.table(filter = "junction", d = SJ_D),
              data.table(filter = "indel", d = INDEL_D))
keep_of <- function(x, f, dd) if (f == "none") rep(TRUE, nrow(x)) else if (f == "junction") is.na(x$d_sj) | x$d_sj > dd else
  is.na(x$d_indel) | x$d_indel > dd
fc <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  f <- grid$filter[i]; dd <- grid$d[i]
  a <- calls[keep_of(calls, f, dd), .(n_calls = .N, n_ag = sum(ag)), by = config]
  s <- f1_of(ct[keep_of(ct, f, dd)])
  s[, `:=`(filter = f, d = dd)]
  merge(s, a, by = "config")
}))
fc <- merge(ranked, fc, by = "config")
fwrite(fc, file.path(OUT, "filter_compare.tsv"), sep = "\t")
b0 <- fc[filter == "none", .(config, gate, t, n_ag0 = n_ag, f1_0 = f1)]
x <- merge(fc[filter != "none"], b0, by = c("config", "gate", "t"))
x[, `:=`(ag_removed_pct = 100 * (1 - n_ag / n_ag0), d_f1 = f1 - f1_0)]
summ <- x[, .(ag_removed_min = min(ag_removed_pct), ag_removed_median = median(ag_removed_pct),
              ag_removed_max = max(ag_removed_pct), d_f1_median = median(d_f1), d_f1_min = min(d_f1),
              d_f1_max = max(d_f1), n_raised = sum(d_f1 > 1e-12), n = .N,
              n_configs_raised_any = uniqueN(config[d_f1 > 1e-12])), by = .(filter, d)]
fwrite(summ, file.path(OUT, "filter_summary.tsv"), sep = "\t")
print(summ, digits = 3)
msg("done ->", OUT)
