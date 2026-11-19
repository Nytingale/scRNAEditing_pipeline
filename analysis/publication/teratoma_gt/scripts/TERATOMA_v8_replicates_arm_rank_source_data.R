#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_replicates_arm_rank_source_data.R   (2026-10-15; split from the _plot.R script 2026-10-25)
# Source data for response figure 17: replicate-aware vs pooled JACUSA2 at the dataset (arm) tier, recall over
# positives with >= 10 reads in at least 2 of T1-T4 NTC.
#   a  replicate-aware vs pooled JACUSA2 F1, the 2 configurations on the REPL_BAMS BAM type(s) x 9 truths
#   b  dataset-tier F1 of the ranked configurations and of replicate-aware JACUSA2 (REPL_BAMS; default
#      sampleUmiMd = UMI-aware MarkDuplicates)
# The variants for positives covered in all 4 teratomas (REPL_SET=reach / all) are in
# ../_superseded/pre_split_2026-10-25/scripts/TERATOMA_v8_replicates_arm_rank_plot.R.
# Input : human/pm/gt_eval_v8/replicates/arm_rank/arm_rank_points.tsv (human/scripts/replicates_arm_rank_v8.R)
# Output: source_data/pm_gt/TERATOMA_v8_replicates_arm_rank.{a_f1_pairs,b_pooled,b_replicate}.tsv
#         (drawn by TERATOMA_v8_replicates_arm_rank_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_replicates_arm_rank"

# REPL_ARM_RANK: 76 configurations, v9 truths = human/pm/gt_eval_v9/replicates/arm_rank/arm_rank_points.tsv
# (human/scripts/replicates_reach_v9.R; with TERATOMA_DATASET=pm_gt76)
SRC <- Sys.getenv("REPL_ARM_RANK", file.path(REPO_DIR, "human", "pm", "gt_eval_v8", "replicates", "arm_rank", "arm_rank_points.tsv"))
stopifnot(file.exists(SRC))
SET <- "reach2"
TRUTH_ORDER <- c(paste("REDIportal", c("cov10_alt1", "alt3", "alt3_e05")),
                 paste("KO", rep(c("depth", "expected"), each = 3), c("t25", "t50", "t75")))
FAM_LAB <- c(REDIportal = "REDIportal truths (3)", KO = "KO truths (6)")

d <- fread(SRC)[set == SET]
stopifnot(nrow(d) > 0)
d[, `:=`(Truth = factor(truth, TRUTH_ORDER), method = caller_factor(caller),
         Fam = factor(FAM_LAB[fam], FAM_LAB))]
TRUTH_KEY <- c(paste("KO expected", c("t25", "t50", "t75")), paste("KO depth", c("t25", "t50", "t75")),
               paste("REDIportal", c("cov10_alt1", "alt3", "alt3_e05")))
TRUTH_SHORT <- setNames(c(sprintf("KO, \u2265 %d%%, power-gated", c(25, 50, 75)),
                          sprintf("KO, \u2265 %d%%, depth-gated", c(25, 50, 75)),
                          "REDIportal \u2265 1 variant read", "REDIportal \u2265 3 variant reads",
                          "REDIportal \u2265 3 variant reads, \u2265 5% editing level"), TRUTH_KEY)
stopifnot(setequal(d$truth, TRUTH_KEY))
d[, TruthY := factor(TRUTH_SHORT[truth], rev(TRUTH_SHORT))]
pool <- d[strategy == "pooled"]
# replicate-aware calls on the BAM types in REPL_BAMS (comma-separated; the response shows UMI-aware MarkDuplicates only)
DEDUP_LAB <- c(sampleUmiDup = "umi_tools", sampleUmiMd = "UMI-aware MarkDuplicates")
REPL_BAMS <- strsplit(Sys.getenv("REPL_BAMS", "sampleUmiMd"), ",", fixed = TRUE)[[1]]
stopifnot(all(REPL_BAMS %in% names(DEDUP_LAB)))
rep <- d[strategy == "jacusa2rep" & bam_type %in% REPL_BAMS]
stopifnot(nrow(rep) > 0)
rep[, Dedup := factor(DEDUP_LAB[bam_type], DEDUP_LAB[names(DEDUP_LAB) %in% bam_type])]

sc <- merge(rep[, .(Fam, truth, config, Dedup, jacusa2rep = f1)],
            pool[, .(truth, config, pooled = f1)], by = c("truth", "config"))
print(sc[, .(n = .N, rep_below_pooled = sum(jacusa2rep < pooled), median_diff = median(jacusa2rep - pooled),
             diff = paste(sprintf("%.3f", range(jacusa2rep - pooled)), collapse = " to ")), by = .(Dedup, Fam)])
print(rep[, .(best_rank = min(rank), worst_rank = max(rank), n_first = sum(rank == 1), .N), by = .(Dedup, Fam)])

write_sd(sc, SRC_DIR, STEM, "a_f1_pairs")
write_sd(pool, SRC_DIR, STEM, "b_pooled")
write_sd(rep, SRC_DIR, STEM, "b_replicate")
