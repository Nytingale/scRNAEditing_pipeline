#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_replicates_source_data.R   (2026-10-15; split from the _plot.R script 2026-10-25)
# Replicate-aware calling vs pooled calling (review R3): the NTC arm called once on the pooled BAM, per teratoma
# then intersected (sites called in >= 2, >= 3 or all 4 of T1-T4), or by JACUSA2 call-1 with the four teratoma
# BAMs as replicates. Scored against the REDIportal truths and the KO truths.
#   f1_pairs  JACUSA2 replicate-aware vs pooled F1, one row per configuration x truth, deduplicated BAMs; recall over
#             all positives (pooled-arm coverage) and over positives with >= 10 reads in each of T1-T4
#   Printed only (numbers in the R3 text): per-strategy medians of the per-configuration means over the truths of a
#   family, deduplicated BAMs (sampleRmDup, sampleUmiDup), all callers.
# Input : human/pm/gt_eval_v8/replicates/summary/replicate_scores.tsv (human/scripts/replicates_summary_v8.R),
#         human/pm/gt_eval_v8/replicates/reachable/reachable_jacusa2.tsv (human/scripts/replicates_reachable_v8.R)
# Output: source_data/pm_gt/TERATOMA_v8_replicates.f1_pairs.tsv (drawn by TERATOMA_v8_replicates_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_replicates"

SRC <- file.path(REPO_DIR, "human", "pm", "gt_eval_v8", "replicates", "summary", "replicate_scores.tsv")
SRC_R <- file.path(REPO_DIR, "human", "pm", "gt_eval_v8", "replicates", "reachable", "reachable_jacusa2.tsv")
stopifnot(file.exists(SRC), file.exists(SRC_R))
DEDUP <- c("sampleRmDup", "sampleUmiDup")
STRAT_LAB <- c(pooled = "Pooled", k2 = "\u22652 teratomas", k3 = "\u22653 teratomas", k4 = "All 4\nteratomas",
               jacusa2rep = "JACUSA2\nreplicate-aware")
FAM_LAB <- c(REDI = "REDIportal truths (3)", KO = "KO truths (6)")

d <- fread(SRC)
d[, fam := factor(FAM_LAB[fifelse(grepl("^REDI", truth), "REDI", "KO")], FAM_LAB)]

a <- d[bam_type %in% DEDUP, .(Precision = mean(precision, na.rm = TRUE), Recall = mean(recall), F1 = mean(f1)),
       by = .(fam, strategy, config, caller)]
a <- melt(a, id.vars = c("fam", "strategy", "config", "caller"), variable.name = "Metric")
a[, `:=`(Strategy = factor(STRAT_LAB[strategy], STRAT_LAB), method = caller_factor(caller))]
print(a[, .(median = round(median(value, na.rm = TRUE), 3), n = .N), by = .(fam, Metric, Strategy)][order(fam, Metric, Strategy)])

SET_LAB <- c(all = "All positives",
             reach = "Positives with \u2265 10 reads\nin every teratoma")
b <- dcast(d[caller == "jacusa2" & strategy %in% c("pooled", "jacusa2rep")],
           fam + truth + config + bam_type ~ strategy, value.var = "f1")[, set := "all"]
r <- fread(SRC_R)[, .(truth, config, bam_type, pooled = f1_pooled, jacusa2rep = f1_jacusa2rep, set = "reach")]
r[, fam := factor(FAM_LAB[fifelse(grepl("^REDI", truth), "REDI", "KO")], FAM_LAB)]
b <- rbind(b, r, use.names = TRUE)[bam_type %in% DEDUP]
PROC_LEVELS <- unname(label_bamType[DEDUP])
b[, `:=`(proc = factor(label_bamType[bam_type], PROC_LEVELS), Set = factor(SET_LAB[set], SET_LAB))]
print(b[, .(n = .N, below_diag = sum(jacusa2rep < pooled)), by = .(Set, fam, proc)])
write_sd(b, SRC_DIR, STEM, "f1_pairs")
