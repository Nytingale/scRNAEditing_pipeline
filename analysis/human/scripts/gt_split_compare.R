#!/usr/bin/env Rscript
# =============================================================================
# gt_split_compare.R   (2026-10-13)
# Split-by-cells test: configuration ranks (gt_rank_v7.R rule: F1 rank per tier, ties min; rank sum over
# arm / sample / celltype; overall rank by (rank sum, -mean F1)) for three scorings, on the configurations
# of the main v8 ranking (mismatch state; the split run adds cellsnp configurations v6 did not report):
#   main       all NTC cells: calls and truth from the same reads (gt_eval_v8/scoring/mismatch)
#   B_full     half-B calls vs the full v8 truth (truth still contains the half-B reads)
#   B_A        half-B calls vs the half-A truth  (truth and calls from disjoint cells)
# B_full vs B_A isolates the shared-reads effect at equal depth; main vs B_full is the depth effect.
# Output: gt_eval_v8/split_cells/
#   compare_config_rank.tsv   gate, t, config, caller, ..., rank / rank_sum / mean_f1 / f1_<tier> per scoring
#   compare_summary.tsv       gate, t, pair, spearman, n_config, top config per scoring, top-5 overlap
#   truth_overlap.tsv         gate, t: positives in A / full / both, informative negatives
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
V8 <- file.path(PM, "gt_eval_v8"); S <- file.path(V8, "split_cells")
TIERS <- c("arm", "sample", "celltype"); KEY <- c("gate", "t")
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
SCN <- Sys.getenv("SCN", "scoring"); SFX <- sub("^scoring", "", SCN)   # SCN=scoring_pwcov: headline denominator
SC <- c(main = file.path(V8, SCN, "mismatch"), B_full = file.path(S, SCN, "truthFull"),
        B_A = file.path(S, SCN, "truthA"))

rd <- function(d) rbindlist(lapply(TIERS, function(tr) fread(file.path(d, tr, "gt_v7_scores_reported.tsv"))[, tier := tr]), fill = TRUE)
sc <- lapply(SC, rd)
cfg <- Reduce(intersect, lapply(sc, function(s) unique(s$config)))
for (nm in names(sc)) message(nm, ": ", uniqueN(sc[[nm]]$config), " configurations, ", length(cfg), " shared")
stopifnot(length(cfg) >= 0.95 * uniqueN(sc$main$config))

rank_of <- function(s) {
  s <- s[config %chin% cfg, c(KEY, "config", COMBO, "tier", "f1"), with = FALSE]
  s[, f1_rank := frank(-fcoalesce(f1, -Inf), ties.method = "min"), by = c(KEY, "tier")]
  w <- dcast(s, as.formula(paste(paste(c(KEY, "config"), collapse = " + "), "~ tier")), value.var = "f1")
  setnames(w, TIERS, paste0("f1_", TIERS))
  r <- s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = c(KEY, "config", COMBO)]
  r[, rank := frank(list(-mean_f1, rank_sum), ties.method = "min"), by = KEY]
  merge(r, w, by = c(KEY, "config"))
}
rk <- lapply(sc, rank_of)
m <- Reduce(function(a, b) merge(a, b, by = c(KEY, "config", COMBO)),
            lapply(names(rk), function(nm) { x <- copy(rk[[nm]]); v <- setdiff(names(x), c(KEY, "config", COMBO))
              setnames(x, v, paste0(v, "_", nm)); x }))
setorderv(m, c(KEY, "rank_main"))
fwrite(m, file.path(S, paste0("compare_config_rank", SFX, ".tsv")), sep = "\t")

PAIRS <- list(c("main", "B_full"), c("main", "B_A"), c("B_full", "B_A"))
top <- function(x, col, k = 1) x[order(get(col))][seq_len(k), config]
sm <- m[, rbindlist(lapply(PAIRS, function(p) {
  a <- paste0("rank_", p[1]); b <- paste0("rank_", p[2])
  data.table(pair = paste(p, collapse = "_vs_"), spearman = cor(.SD[[a]], .SD[[b]], method = "spearman"), n_config = .N,
             top_x = top(.SD, a), top_y = top(.SD, b),
             top5_shared = length(intersect(top(.SD, a, 5), top(.SD, b, 5))))
})), by = KEY]
fwrite(sm, file.path(S, paste0("compare_summary", SFX, ".tsv")), sep = "\t")
print(sm[, .(gate, t, pair, spearman = round(spearman, 3), top_x, top_y, top5_shared)])

ta <- fread(file.path(S, "truth", "v7_nomm_pgp1", "gt_v7_truth_sites.tsv.gz"))
tf <- fread(file.path(V8, "truth", "gt_v7_truth_sites.tsv.gz"))
ov <- rbindlist(lapply(c("t25", "t50", "t75"), function(tn) {
  cc <- paste0("class_", tn)
  merge(ta[, .(pos_A = sum(get(cc) == "positive"), neg_A = sum(get(cc) == "informative_negative")), by = gate],
        tf[, .(pos_full = sum(get(cc) == "positive"), neg_full = sum(get(cc) == "informative_negative")), by = gate], by = "gate")[
    , `:=`(t = tn, pos_both = sapply(gate, function(g) length(intersect(ta[gate == g & get(cc) == "positive", site],
                                                                         tf[gate == g & get(cc) == "positive", site]))))]
}))
setcolorder(ov, c("gate", "t"))
fwrite(ov, file.path(S, "truth_overlap.tsv"), sep = "\t")
print(ov)
