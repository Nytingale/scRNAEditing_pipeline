#!/usr/bin/env Rscript
# replicates_ranksum_v8.R   (2026-10-15; Reviewer 3 replicates)
# What-if rank sum for replicate-aware JACUSA2: each of its 4 deduplicated-BAM configurations is added as a 59th
# configuration whose dataset-tier F1 is its own and whose sample- and cell-type-tier F1 are those of the same pooled
# JACUSA2 configuration (replicate mode has no sample/cell-type equivalent). F1 ranks within tier (ties min), rank sum,
# overall rank by (rank sum, -mean F1), as gt_rank_v7.R. Dataset-tier recall: "all" = all positives (main ranking);
# "reach" / "reach2" = positives with >= 10 reads in each of / at least 2 of T1-T4, for all 59 configurations
# (sample/cell-type tiers unchanged).
# Output: gt_eval_v8/replicates/arm_rank/ranksum_jacusa2rep.tsv
suppressPackageStartupMessages(library(data.table))
G <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
TIERS <- c("arm", "sample", "celltype")

ko <- rbindlist(lapply(TIERS, function(tr)
  fread(file.path(G, "scoring", "mismatch", tr, "gt_v7_scores_reported.tsv"))[, .(truth = paste("KO", gate, t), config, tier = tr, f1)]))
redi <- rbindlist(lapply(c("cov10_alt1", "alt3", "alt3_e05"), function(tt)
  fread(file.path(G, "rediportal", "scoring_mismatch", paste0("gt_tier_scoring_", tt, ".tsv")))[
    , .(truth = paste("REDIportal", tt), config = paste(caller, bam_type, strand_mode, snp_subtract, sep = "__"),
        tier = c(dataset = "arm", sample = "sample", celltype = "celltype")[tier], f1)]))
s <- rbind(ko, redi)
stopifnot(s[, .N, by = .(truth, tier)]$N == 58)

pts <- fread(file.path(G, "replicates", "arm_rank", "arm_rank_points.tsv"))
chk <- merge(pts[set == "all" & strategy == "pooled", .(truth, config, f1)], s[tier == "arm"], by = c("truth", "config"))
stopifnot(nrow(chk) == 58 * 9, max(abs(chk$f1.x - chk$f1.y)) < 1e-9)
reps <- pts[strategy == "jacusa2rep", unique(config)]

out <- rbindlist(lapply(unique(pts$set), function(st) rbindlist(lapply(reps, function(r) {
  base <- copy(s)
  if (st != "all") {
    base <- merge(base, pts[set == st & strategy == "pooled", .(truth, config, f1r = f1)], by = c("truth", "config"),
                  all.x = TRUE)
    base[tier == "arm", f1 := f1r][, f1r := NULL]
  }
  add <- base[config == r][, strategy := "jacusa2rep"]
  add[tier == "arm", f1 := pts[set == st & strategy == "jacusa2rep" & config == r][match(add[tier == "arm", truth], truth), f1]]
  x <- rbind(base[, strategy := "pooled"], add)
  x[, f1_rank := frank(-f1, ties.method = "min"), by = .(truth, tier)]
  rs <- x[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1)), by = .(truth, strategy, config)]
  rs[, overall_rank := frank(list(-mean_f1, rank_sum), ties.method = "min"), by = truth]
  w <- dcast(x, truth + strategy + config ~ tier, value.var = "f1_rank")
  rs <- merge(rs, w, by = c("truth", "strategy", "config"))
  top <- rs[overall_rank == 1, .(truth, top_config = config, top_rank_sum = rank_sum)]
  pj <- rs[strategy == "pooled" & config == r, .(truth, pooled_same_overall = overall_rank)]
  Reduce(function(a, b) merge(a, b, by = "truth"),
         list(rs[strategy == "jacusa2rep"], pj, top))[, set := st]
}))))
setcolorder(out, c("set", "truth", "config", "overall_rank", "rank_sum", "arm", "sample", "celltype", "mean_f1",
                   "pooled_same_overall", "top_config", "top_rank_sum"))
fwrite(out[, !"strategy"], file.path(G, "replicates", "arm_rank", "ranksum_jacusa2rep.tsv"), sep = "\t")
print(out[, .(best = min(overall_rank), worst = max(overall_rank), n_first = sum(overall_rank == 1), .N), by = set])
print(out[, .(set, truth, config = sub("jacusa2__", "", config), overall_rank, rank_sum, arm, sample, celltype,
              pooled_same_overall, top_rank_sum)][order(set, truth, overall_rank)], nrows = 100)
