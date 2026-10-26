#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_negatives_source_data.R   (2026-10-15; split from TERATOMA_v8_negatives_plot.R 2026-10-25)
# Source data for response figure 3: pooled NTC editing level of the informative negatives vs the positives
# (t = 50%), per coverage gate, and the share of each class in REDIportal.
# Input : human/pm/gt_eval_v8/truth/negatives_describe_sites.tsv.gz (human/scripts/gt_v8_negatives_describe.R)
# Output: source_data/pm_gt/TERATOMA_v8_negatives_ntc_editing.{sites,labels}.tsv
#         (drawn by TERATOMA_v8_negatives_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_negatives_ntc_editing"

SRC <- file.path(V8_DIR, "truth", "negatives_describe_sites.tsv.gz")
stopifnot(file.exists(SRC))
d <- fread(SRC)[class %in% c("informative_negative", "positive_t50")]
CLASS_LAB <- c(positive_t50 = "Positives (KO depletion \u2265 50%)", informative_negative = "Informative negatives")
d[, `:=`(Class = factor(CLASS_LAB[class], CLASS_LAB), Gate = factor(V8_GATE_LAB[gate], V8_GATE_LAB))]
lab <- d[, .(n = .N, med = median(edit), redi = mean(in_redi)), by = .(Gate, Class)]
lab[, txt := sprintf("n = %s, median %.2f, %.0f%% are REDIportal sites", format(n, big.mark = ",", trim = TRUE), med, 100 * redi)]

write_sd(d[, .(Gate, Class, site, edit, in_redi)], SRC_DIR, STEM, "sites")
write_sd(lab, SRC_DIR, STEM, "labels")
