#!/usr/bin/env Rscript
# gt_v9_unit_cov_collect.R   (2026-11-07)
# Per-unit coverage on the read strand of each v9 positive (gt_v9_build.sh collect): forward pileup
# (gt_v8_unit_cov.sh LIB=forward) for read A>G sites, reverse pileup for read T>C sites.
# -> GT_V9_DIR/unit_cov/ko_truth_unit_cov.tsv.gz (unit site coverage variant_cov; gt_score_v7.R GT_UNIT_COV)
suppressPackageStartupMessages(library(data.table))
PM = file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
V9 = Sys.getenv("GT_V9_DIR", file.path(PM, "gt_eval_v9"))
tr = fread(file.path(V9, "truth/gt_v9_positive_sites.tsv.gz"))
rt = unique(tr[, .(site, read_type)])
stopifnot(rt[, .N, by = site][, all(N == 1L)], all(rt$read_type %in% c("AG", "TC")))
f = fread(file.path(V9, "unit_cov_fwd/ko_truth_unit_cov.tsv.gz"))
r = fread(file.path(V9, "unit_cov_rev/ko_truth_unit_cov.tsv.gz"))
uc = rbind(f[site %chin% rt[read_type == "AG", site]], r[site %chin% rt[read_type == "TC", site]])
dir.create(file.path(V9, "unit_cov"), showWarnings = FALSE)
fwrite(uc, file.path(V9, "unit_cov/ko_truth_unit_cov.tsv.gz"), sep = "\t")
if (Sys.getenv("GT_V9_UCOV_PARITY", "1") == "1") {   # 0 for BAMs other than the v8 ones (seed replicates)
old = fread(file.path(PM, "gt_eval_v8/unit_cov/ko_truth_unit_cov.tsv.gz"))
m = merge(old, f, by = c("unit", "site"))
message("forward pileup vs v8 unit_cov at shared sites: ", nrow(m), " rows, identical coverage ",
        sum(m$coverage.x == m$coverage.y), ", variant_cov ", sum(m$variant_cov.x == m$variant_cov.y))
stopifnot(all(m$coverage.x == m$coverage.y), all(m$variant_cov.x == m$variant_cov.y))
}
print(uc[unit == "NTC_ctrl", .(sites = .N, cov10 = sum(coverage >= 10)), by = .(read_TC = site %chin% rt[read_type == "TC", site])])
