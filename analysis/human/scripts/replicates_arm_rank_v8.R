#!/usr/bin/env Rscript
# replicates_arm_rank_v8.R   (2026-10-15; Reviewer 3 replicates)
# Replicate-aware JACUSA2 is not a ranked configuration (the ranking sums arm, sample and cell-type ranks, and only the
# arm has replicates). Here its 4 deduplicated-BAM configurations are placed among the 58 ranked configurations at the
# dataset (arm) tier: per truth, F1 rank = 1 + number of ranked configurations with higher F1. Recall over all positives
# (replicate_scores.tsv), over positives with >= 10 reads in each of T1-T4 (reachable/reachable_scores.tsv; "reach") and,
# when present, in at least 2 of T1-T4 (reachable_min2/; "reach2").
# Output: gt_eval_v8/replicates/arm_rank/{arm_rank_points.tsv, arm_rank_jacusa2rep.tsv, arm_rank_top.tsv}
suppressPackageStartupMessages(library(data.table))
R <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8/replicates")
O <- file.path(R, "arm_rank"); dir.create(O, showWarnings = FALSE)
DEDUP <- c("sampleRmDup", "sampleUmiDup")

all <- fread(file.path(R, "summary", "replicate_scores.tsv"))[strategy %in% c("pooled", "jacusa2rep")]
ranked <- unique(all[strategy == "pooled", config])
stopifnot(length(ranked) == 58)
REACH <- c(reach = "reachable", reach2 = "reachable_min2")
REACH <- REACH[file.exists(file.path(R, REACH, "reachable_scores.tsv"))]
pts <- rbind(all[, .(set = "all", strategy, config, caller, bam_type, truth, f1)],
             rbindlist(lapply(names(REACH), function(st)
               fread(file.path(R, REACH[[st]], "reachable_scores.tsv"))[strategy %in% c("pooled", "jacusa2rep"),
                 .(set = st, strategy, config, caller, bam_type, truth, f1)])))
pts <- pts[(strategy == "pooled" & config %in% ranked) | (strategy == "jacusa2rep" & bam_type %in% DEDUP)]
stopifnot(pts[strategy == "pooled", .N, by = .(set, truth)]$N == 58)
pts[, fam := fifelse(grepl("^REDI", truth), "REDIportal", "KO")]

pool <- pts[strategy == "pooled"]
rep <- pts[strategy == "jacusa2rep"]
rep[, rank := mapply(function(s, t, f) 1L + pool[set == s & truth == t, sum(f1 > f)], set, truth, f1)]
pool[, rank := frank(-f1, ties.method = "min"), by = .(set, truth)]
top <- pool[rank == 1, .(top_config = paste(config, collapse = ";"), top_f1 = f1[1]), by = .(set, fam, truth)]
rep <- merge(rep, top[, .(set, truth, top_f1)], by = c("set", "truth"))
rep[, f1_gap_to_top := f1 - top_f1]
fwrite(rbind(pool, rep[, names(pool), with = FALSE]), file.path(O, "arm_rank_points.tsv"), sep = "\t")
fwrite(rep, file.path(O, "arm_rank_jacusa2rep.tsv"), sep = "\t")
fwrite(top, file.path(O, "arm_rank_top.tsv"), sep = "\t")
print(rep[, .(best_rank = min(rank), worst_rank = max(rank), rank1 = sum(rank == 1), n = .N,
              gap = paste(sprintf("%.3f", range(f1_gap_to_top)), collapse = " to ")), by = .(set, fam)])
print(top)
