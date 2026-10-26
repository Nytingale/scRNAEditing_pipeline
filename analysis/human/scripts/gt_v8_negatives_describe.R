#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_negatives_describe.R   (2026-10-15)
#
# Describes the v8 informative negatives against the positives (review response, Reviewer 3 ASE point):
# pooled NTC editing level (NTC_ctrl unit, UMI-dedup stranded pileup at the truth positions) and the
# fraction in REDIportal.
# Inputs : gt_eval_v8/truth/gt_v7_truth_sites.tsv.gz, gt_eval_v8/unit_cov/ko_truth_unit_cov.tsv.gz,
#          gt_eval_v6/rediportal/bed/chr*.bed (0-based start, end = 1-based position)
# Output : gt_eval_v8/truth/negatives_describe.tsv (gate, class, n, median_edit, frac_edit_ge_0.9, frac_rediportal)
#          gt_eval_v8/truth/negatives_describe_sites.tsv.gz (gate, site, class, edit, in_redi; for TERATOMA_v8_negatives_plot.R)
# Usage  : Rscript gt_v8_negatives_describe.R   (login node OK)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM = file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
V8 = file.path(PM, "gt_eval_v8")

tr = fread(file.path(V8, "truth/gt_v7_truth_sites.tsv.gz"))
cls = rbind(
  tr[class_t25 == "informative_negative", .(gate, site, chr, pos, class = "informative_negative")],
  rbindlist(lapply(c("t25", "t50", "t75"), function(t_)
    tr[get(paste0("class_", t_)) == "positive", .(gate, site, chr, pos, class = paste0("positive_", t_))])))
stopifnot(identical(tr[class_t25 == "informative_negative", sort(site)], tr[class_t75 == "informative_negative", sort(site)]))

cov = fread(file.path(V8, "unit_cov/ko_truth_unit_cov.tsv.gz"))[unit == "NTC_ctrl", .(site, edit = variant_cov / coverage)]
cls = merge(cls, cov, by = "site", all.x = TRUE)
stopifnot(!anyNA(cls$edit))

redi = rbindlist(lapply(list.files(file.path(PM, "gt_eval_v6/rediportal/bed"), "^chr.*\\.bed$", full.names = TRUE), function(f)
  fread(f, header = FALSE, select = c(1, 3), col.names = c("chr", "pos"))[pos %in% cls$pos]))
redi = unique(redi)[, in_redi := TRUE]
cls = merge(cls, redi, by = c("chr", "pos"), all.x = TRUE)[is.na(in_redi), in_redi := FALSE]

out = cls[, .(n = .N, median_edit = round(median(edit), 3), frac_edit_ge_0.9 = round(mean(edit >= 0.9), 3),
              frac_rediportal = round(mean(in_redi), 3)), by = .(gate, class)][order(gate, class)]
fwrite(out, file.path(V8, "truth/negatives_describe.tsv"), sep = "\t")
fwrite(cls[, .(gate, site, class, edit, in_redi)], file.path(V8, "truth/negatives_describe_sites.tsv.gz"), sep = "\t")
print(out)
