#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_sj_indel_calls.R   (2026-10-23; Reviewer 2, major 3: false-positive SNVs near splice junctions and indels)
# Are the called sites near splice junctions or indels enriched for false positives? Merged NTC arm tier, the 58
# ranked configurations, main call state (step 5: after removal of sites near other mismatch types), all 12 types.
#   d_sj     distance to the nearest splice junction boundary (first or last intron base of a STAR junction with >= 5
#            uniquely mapped reads summed over the 8 teratoma runs, or annotated), as gt_v8_sj_homopolymer_truth.R
#   d_indel  distance to the nearest indel in the merged NTC reads with >= MIN_INDEL reads (ntc_indel_positions_v8.sh);
#            deletion = its first and last deleted base, insertion = the two reference bases around it
# False-positive measures for every call:
#   mismatch type   A-to-G (RNA frame) vs the other 11 types, which ADAR does not make
#   KO response     evaluable (NTC cov >= 10, m = KO cov x NTC variant fraction >= 3 in either KO) and responding
#                   (KO variant fraction <= 0.5 x NTC) as gt_v8_callset_ko_response.R
#   truth           TP / FP at the KO truth sites as gt_score_v7.R (TP: correct call at a positive; FP: wrong-type call
#                   at a positive, or any call at an informative negative)
# Filters: remove calls with d_sj <= d (d = 1-6, 8, 10) or d_indel <= d (d = 1, 2, 3, 5, 10); per configuration the
# A-to-G sites, A-to-G %, KO response and arm F1 (headline recall denominator, GT_UNIT_RULE=gate_cov) are recomputed.
# Outputs (human/pm/gt_eval_v8/indel_sj/):
#   call_sites.tsv.gz    unique called sites with d_sj, d_indel, KO response, truth labels
#   distance_bins.tsv    configuration x feature x distance bin: calls, A-to-G, KO response, TP, FP
#   filter_compare.tsv   configuration x filter x d (d = 0: no filter): sites, A-to-G %, KO response, F1 per truth
# =============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("SCRNA_DIR", ".")
PM <- file.path(ROOT, "human/pm")
V8 <- file.path(PM, "gt_eval_v8")
SP <- file.path(PM, "gt_eval_v6/spectrum12")
RES <- file.path(PM, "results_gt_v8/NTC_ctrl/04_filtered")
OUT <- file.path(V8, "indel_sj")
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
MIN_INDEL <- as.integer(Sys.getenv("MIN_INDEL", "3"))
M_EXP <- 3; MIN_COV <- 10; T_DROP <- 0.5
SJ_D <- c(1:6, 8L, 10L)
INDEL_D <- c(1L, 2L, 3L, 5L, 10L)
BRK <- c(-1, 5, 10, 20, 50, 100, 500, Inf)
LAB <- c("0-5", "6-10", "11-20", "21-50", "51-100", "101-500", "> 500")
COMP <- c(A = "T", C = "G", G = "C", T = "A")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

ranked <- unique(fread(file.path(V8, "state_summary/config_ranks.tsv"))[, c(COMBO, "config"), with = FALSE])
stopifnot(nrow(ranked) == 58L)
fs <- list.files(RES, recursive = TRUE, full.names = TRUE, pattern = "_step5_nomm\\.tsv\\.gz$")
meta <- data.table(f = fs)
meta[, c("caller", "bam_type", "strand_mode", "snp_dir") := tstrsplit(sub(paste0(RES, "/"), "", f, fixed = TRUE), "/",
                                                                         fixed = TRUE, keep = 1:4)]
meta[, snp_subtract := sub("snp_", "", snp_dir)]
meta[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
meta <- meta[config %in% ranked$config]
stopifnot(nrow(meta) == 58L)
calls <- unique(rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  x <- fread(meta$f[i], select = c("chr", "pos", "varType", "strand"), showProgress = FALSE)
  x[, config := meta$config[i]]
})))
calls[, site := paste0(chr, ":", pos)]
msg("calls:", nrow(calls))
sites <- unique(calls[, .(chr, pos, site, strand, varType)])

# ---- distances ------------------------------------------------------------------
nearest <- function(pts, q) {  # pts: chr, b; q: chr, pos -> distance to nearest b on the chromosome
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
upos <- unique(sites[, .(chr, pos)])
upos[, d_sj := nearest(bd, upos)]

inf <- list.files(file.path(OUT, "indels"), "\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(inf) == 23)
ind <- rbindlist(lapply(inf, fread, col.names = c("chr", "pos", "type", "len", "n_reads")))
ind <- ind[, .(n_reads = sum(n_reads)), by = .(chr, pos, type, len)][n_reads >= MIN_INDEL]
msg("indel position x length rows with >=", MIN_INDEL, "reads:", nrow(ind), "(", ind[type == "D", .N], "deletions)")
ib <- unique(rbind(ind[type == "D", .(chr, b = pos)], ind[type == "D", .(chr, b = pos + len - 1L)],
                   ind[type == "I", .(chr, b = pos - 1L)], ind[type == "I", .(chr, b = pos)]))
upos[, d_indel := nearest(ib, upos)]
sites <- upos[sites, on = .(chr, pos)]

# ---- KO response (gt_v8_callset_ko_response.R) --------------------------------------
pidx <- fread(file.path(SP, "positions_noSNP.bed"), header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
need <- unique(sites[, .(chr, pos)])[pidx, on = .(chr, pos), nomatch = NULL]
need[, idx := .I]
cols <- c("ref_plus", "ref_minus", "A", "C", "G", "T", "a", "c", "g", "t")
read_unit <- function(u) {
  f <- list.files(file.path(SP, "pileups"), pattern = paste0("^", u, "__.*\\.tsv\\.gz$"), full.names = TRUE)
  d <- rbindlist(lapply(f, function(x) fread(cmd = paste("gzip -cd", shQuote(x)), sep = "\t", showProgress = FALSE)))
  d <- need[d, on = .(chr, pos), nomatch = NULL]
  m <- matrix(0L, nrow(need), length(cols), dimnames = list(NULL, cols))
  for (cn in cols) m[d$idx, cn] <- as.integer(d[[cn]])
  msg(u, ":", nrow(d), "positions")
  m
}
U <- lapply(c(ntc = "merged_NTC_ctrl", ko1 = "merged_ADAR1_KO", ko2 = "merged_ADARB1_KO"), read_unit)
sites <- need[sites, on = .(chr, pos)]
sites[, `:=`(Y = substr(varType, 2, 2), plus = strand == "+")]
ok <- !is.na(sites$idx)
counts <- function(m) {
  alt <- refn <- rep(NA_integer_, nrow(sites))
  i <- sites$idx[ok]; pl <- sites$plus[ok]; Y <- sites$Y[ok]
  acol <- ifelse(pl, match(Y, cols), match(tolower(COMP[Y]), cols))
  alt[ok] <- m[cbind(i, acol)]
  refn[ok] <- ifelse(pl, m[i, "ref_plus"], m[i, "ref_minus"])
  list(alt = alt, cov = refn + alt)
}
S <- lapply(U, counts)
rm(U); gc()
fr <- function(s) fifelse(s$cov > 0, s$alt / pmax(s$cov, 1), NA_real_)
sites[, `:=`(ntc_cov = S$ntc$cov, ntc_alt = S$ntc$alt, ntc_vf = fr(S$ntc),
             ko1_cov = S$ko1$cov, ko1_vf = fr(S$ko1), ko2_cov = S$ko2$cov, ko2_vf = fr(S$ko2))]
ev <- function(kc) !is.na(kc) & sites$ntc_cov >= MIN_COV & sites$ntc_alt > 0 & kc * sites$ntc_vf >= M_EXP
sites[, `:=`(ev1 = ev(ko1_cov), ev2 = ev(ko2_cov))]
sites[, `:=`(evaluable = ev1 | ev2,
             responds = (ev1 & ko1_vf <= (1 - T_DROP) * ntc_vf) | (ev2 & ko2_vf <= (1 - T_DROP) * ntc_vf))]
sites[is.na(responds), responds := FALSE]
sites[, ag := varType == "AG"]

# ---- truth labels (t50, both gates) -------------------------------------------------
truth <- fread(file.path(V8, "truth/gt_v7_truth_sites.tsv.gz"), select = c("gate", "site", paste0("class_t", c(25, 50, 75))))
for (g in c("expected", "depth")) {
  cl <- truth[gate == g, .(site, cls = class_t50)]
  sites[cl, (paste0("cls_", g)) := i.cls, on = "site"]
  sites[, (paste0("tp_", g)) := ag & get(paste0("cls_", g)) %chin% "positive"]
  sites[, (paste0("fp_", g)) := (!ag & get(paste0("cls_", g)) %chin% "positive") |
                                 get(paste0("cls_", g)) %chin% "informative_negative"]
}
fwrite(sites[, .(chr, pos, strand, varType, d_sj, d_indel, ntc_cov, ntc_alt, evaluable, responds, cls_expected,
                 cls_depth)], file.path(OUT, "call_sites.tsv.gz"), sep = "\t")
cs <- sites[, .(chr, pos, strand, varType, d_sj, d_indel, ag, evaluable, responds, tp_expected, fp_expected, tp_depth,
                fp_depth)][calls, on = .(chr, pos, strand, varType)]

# ---- distance bins ------------------------------------------------------------------
bin_of <- function(d) factor(LAB[findInterval(fifelse(is.na(d), Inf, d), BRK, left.open = TRUE)], LAB)
summ <- function(x) x[, .(n_calls = .N, n_ag = sum(ag), n_eval = sum(evaluable), n_resp = sum(responds),
                          n_eval_ag = sum(evaluable & ag), n_resp_ag = sum(responds & ag),
                          n_eval_other = sum(evaluable & !ag), n_resp_other = sum(responds & !ag),
                          tp_expected = sum(tp_expected), fp_expected = sum(fp_expected),
                          tp_depth = sum(tp_depth), fp_depth = sum(fp_depth))]
bins <- rbind(cs[, c(list(feature = "junction"), summ(.SD)), by = .(config, bin = bin_of(d_sj))],
              cs[, c(list(feature = "indel"), summ(.SD)), by = .(config, bin = bin_of(d_indel))])
fwrite(merge(ranked, bins, by = "config"), file.path(OUT, "distance_bins.tsv"), sep = "\t")

# ---- filters: arm F1 recomputed ------------------------------------------------------
uc <- fread(file.path(V8, "unit_cov/ko_truth_unit_cov.tsv.gz"))[unit == "NTC_ctrl"]
uc[, edit := variant_cov / coverage]
covered <- rbind(uc[coverage >= MIN_COV, .(site, gate = "depth")],
                 uc[coverage >= MIN_COV & coverage * edit >= M_EXP, .(site, gate = "expected")])
tl <- melt(truth, id.vars = c("gate", "site"), variable.name = "t", value.name = "cls")[cls != "unresolved"]
tl[, t := sub("class_", "", t)]
tl[, covered := FALSE][covered, covered := TRUE, on = .(site, gate)]
npu <- tl[cls == "positive" & covered == TRUE, .(n_pos_unit = .N), by = .(gate, t)]
ct <- unique(cs[site %chin% truth$site, .(config, site, ag, d_sj, d_indel)])
f1_of <- function(keep_ct) {
  sc <- keep_ct[, .(any_correct = any(ag), any_wrong = any(!ag)), by = .(config, site)]
  x <- merge(sc, tl, by = "site", allow.cartesian = TRUE)
  s <- x[, .(n_tp = sum(cls == "positive" & any_correct), n_tp_cov = sum(cls == "positive" & any_correct & covered),
             n_fp = sum(cls == "positive" & any_wrong) + sum(cls == "informative_negative")), by = .(config, gate, t)]
  s <- merge(CJ(config = ranked$config, gate = c("depth", "expected"), t = c("t25", "t50", "t75")), s,
             by = c("config", "gate", "t"), all.x = TRUE)
  for (k in c("n_tp", "n_tp_cov", "n_fp")) set(s, which(is.na(s[[k]])), k, 0L)
  s <- merge(s, npu, by = c("gate", "t"))
  s[, `:=`(precision = fifelse(n_tp + n_fp > 0, n_tp / (n_tp + n_fp), NA_real_), recall = n_tp_cov / n_pos_unit)]
  s[, f1 := fifelse(is.na(precision) | precision + recall == 0, 0, 2 * precision * recall / (precision + recall))]
  s
}
base_f1 <- f1_of(ct)
ref <- fread(file.path(V8, "scoring_pwcov/mismatch/arm/gt_v7_unit_scores.tsv"))[config %in% ranked$config]
chk <- merge(base_f1, ref[, .(config, gate, t, f1_ref = f1)], by = c("config", "gate", "t"))
stopifnot(nrow(chk) == 58L * 6L)
msg("arm F1 vs published scoring, max |diff|:", signif(chk[, max(abs(f1 - f1_ref))], 3))
stopifnot(chk[, max(abs(f1 - f1_ref))] < 1e-9)

grid <- rbind(data.table(filter = "none", d = 0L), data.table(filter = "junction", d = SJ_D),
              data.table(filter = "indel", d = INDEL_D))
fc <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  f <- grid$filter[i]; dd <- grid$d[i]
  keep <- function(x) if (f == "none") rep(TRUE, nrow(x)) else if (f == "junction") is.na(x$d_sj) | x$d_sj > dd else
    is.na(x$d_indel) | x$d_indel > dd
  k <- cs[keep(cs)]
  a <- k[, .(n_calls = .N, n_ag = sum(ag), ag_pct = 100 * mean(ag), n_eval = sum(evaluable), n_resp = sum(responds),
             n_eval_ag = sum(evaluable & ag), n_resp_ag = sum(responds & ag)), by = config]
  s <- f1_of(ct[keep(ct)])
  w <- dcast(s, config ~ paste(gate, t, sep = "_"), value.var = c("f1", "precision", "recall"))
  merge(a, w, by = "config")[, `:=`(filter = f, d = dd)]
}))
fc <- merge(ranked, fc, by = "config")
fwrite(fc, file.path(OUT, "filter_compare.tsv"), sep = "\t")
base <- fc[filter == "none", .(config, n_ag0 = n_ag, ag_pct0 = ag_pct, f1_0 = f1_expected_t50)]
print(merge(fc, base, by = "config")[, .(sites_removed_pct = median(100 * (1 - n_ag / n_ag0)),
                                          d_ag_pct = median(ag_pct - ag_pct0),
                                          d_f1_expected_t50 = median(f1_expected_t50 - f1_0)), by = .(filter, d)])
