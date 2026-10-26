#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_strong_truth.R   (2026-10-13)
# Sensitivity test (supplement): v8 KO truth restricted to positives with strong NTC evidence, to check
# that the ranking does not hinge on sites near the truth thresholds (where a caller sharing the truth's
# pileup engine, samtools mpileup, would be favoured). Positives with < ALT variant reads (stranded,
# merged NTC UMI-dedup BAM, BQ/MQ >= 30; gt_v8_unit_cov.sh) become "unresolved" (unscored); informative
# negatives are unchanged.
# Usage : Rscript gt_v8_strong_truth.R 10 20
# Output: gt_eval_v8/sensitivity/strong_alt<ALT>/truth/gt_v7_truth_sites.tsv.gz, truth_counts.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
V8 = file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
tr = fread(file.path(V8, "truth", "gt_v7_truth_sites.tsv.gz"))
cv = fread(file.path(V8, "unit_cov", "ko_truth_unit_cov.tsv.gz"))[unit == "NTC_ctrl", .(site, ntc_alt = variant_cov)]
stopifnot(all(tr$site %in% cv$site))
for (alt in as.integer(commandArgs(trailingOnly = TRUE))) {
  x = merge(tr, cv, by = "site", sort = FALSE)
  for (cl in c("class_t25", "class_t50", "class_t75")) x[get(cl) == "positive" & ntc_alt < alt, (cl) := "unresolved"]
  o = file.path(V8, "sensitivity", paste0("strong_alt", alt), "truth"); dir.create(o, recursive = TRUE, showWarnings = FALSE)
  fwrite(x[, names(tr), with = FALSE], file.path(o, "gt_v7_truth_sites.tsv.gz"), sep = "\t")
  n = x[, .(t25 = sum(class_t25 == "positive"), t50 = sum(class_t50 == "positive"), t75 = sum(class_t75 == "positive"),
            neg = sum(class_t50 == "informative_negative")), by = gate][, min_ntc_alt := alt]
  n0 = tr[, .(t25_all = sum(class_t25 == "positive"), t50_all = sum(class_t50 == "positive"), t75_all = sum(class_t75 == "positive")), by = gate]
  n = merge(n, n0, by = "gate"); fwrite(n, file.path(o, "truth_counts.tsv"), sep = "\t"); print(n)
}
