#!/usr/bin/env Rscript
# gt_v9_baq_compare.R   (2026-11-17; response paragraph 382, BCFtools BAQ on vs off)
# Mean F1 (three tiers) of the 20 BCFtools configurations with BAQ on (published, scoring_pwcov76) and off
# (human/pm/scripts/baqoff_v9.sh), under each of the six KO-based v9 truth sets.
# Input : human/pm/gt_eval_v9/{scoring_pwcov76/mismatch,baq_compare/scoring_pwcov}/config_rank.tsv
# Output: human/pm/gt_eval_v9/baq_compare/compare_f1.tsv (config x truth), compare_f1_summary.tsv (per truth)
suppressPackageStartupMessages(library(data.table))
V9 <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v9")
cols <- c("gate", "t", "config", "bam_type", "mean_f1")
on <- fread(file.path(V9, "scoring_pwcov76", "mismatch", "config_rank.tsv"))[caller == "bcftools", ..cols]
off <- fread(file.path(V9, "baq_compare", "scoring_pwcov", "config_rank.tsv"))[, ..cols]
stopifnot(uniqueN(on$config) == 20, setequal(unique(off$config), unique(on$config)))
x <- merge(on, off, by = c("gate", "t", "config", "bam_type"), suffixes = c("_baq_on", "_baq_off"))
stopifnot(nrow(x) == nrow(on))
x[, diff := mean_f1_baq_on - mean_f1_baq_off]
fwrite(x, file.path(V9, "baq_compare", "compare_f1.tsv"), sep = "\t")
s <- x[, .(n_configs = .N, n_raised = sum(diff > 0), median_diff = median(diff), min_diff = min(diff), max_diff = max(diff)),
       by = .(gate, t)]
fwrite(s, file.path(V9, "baq_compare", "compare_f1_summary.tsv"), sep = "\t")
print(s, digits = 3)
cat("\nAll truths: raised in", x[, sum(diff > 0)], "of", nrow(x), "config x truth; configs raised under every truth:",
    x[, all(diff > 0), by = config][, sum(V1)], "of 20; median per truth", round(min(s$median_diff), 3), "to",
    round(max(s$median_diff), 3), "\n")
print(x[, .(median_diff = median(diff)), by = bam_type], digits = 3)
