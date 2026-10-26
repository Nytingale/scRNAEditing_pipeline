#!/usr/bin/env Rscript
# =============================================================================
# replicates_kfilter_v8.R   (2026-10-15; Reviewer 3 replicates; from _v6)
# Replicate filter for every configuration: a site is kept if it is in the main-state call set
# (_step5_nomm) of >= k of the four teratoma NTC samples (T1-T4, results_gt_persample_v8), k = 2, 3, 4.
# varType / strand = the type and strand called in most of those teratomas (ties: first teratoma). Written as arm-level
# call sets in the results-tree layout gt_score_rediportal_v8.R / gt_score_v7.R read (site, varType, strand).
# Output: gt_eval_v8/replicates/results_k<k>_arm/NTC_ctrl/04_filtered/<caller>/<bam>/<strand>/<snp>/
#         NTC_ctrl_<caller>_<bam>_<strand>_snp<x>_step5_nomm.tsv.gz
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
R <- file.path(PM, "gt_eval_v8", "replicates")
S <- file.path(PM, "results_gt_persample_v8")
KS <- 2:4
f <- unlist(lapply(paste0("T", 1:4, "_NTC_ctrl"), function(u)
  list.files(file.path(S, u, "04_filtered"), "_step5_nomm\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)))
fi <- data.table(file = f, unit = sub(".*/(T[1-4]_NTC_ctrl)/04_filtered/.*", "\\1", f),
                 rel = sub(".*/04_filtered/", "", f))
fi[, cfg := dirname(rel)]
fi[, stem := sub("^T[1-4]_NTC_ctrl_", "NTC_ctrl_", basename(rel))]
bad <- fi[, .N, by = cfg][N != 4]
if (nrow(bad)) { print(bad); stop("configurations without all four teratomas") }
cat(uniqueN(fi$cfg), "configurations x 4 teratomas\n")
for (cf in unique(fi$cfg)) {
  x <- rbindlist(lapply(fi[cfg == cf, file], function(p)
    if (file.size(p) == 0) NULL else fread(p, select = c("site", "varType", "strand"))[, f := p]))
  st <- fi[cfg == cf, stem][1]
  if (!nrow(x)) x <- data.table(site = character(), varType = character(), strand = character(), f = character())
  x <- unique(x)
  cnt <- x[, .(n_ter = uniqueN(f)), by = site]
  vt <- x[, .N, by = .(site, varType, strand)][order(site, -N)][, .SD[1], by = site][, .(site, varType, strand)]
  for (k in KS) {
    out <- file.path(R, paste0("results_k", k, "_arm"), "NTC_ctrl", "04_filtered", cf)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    fwrite(merge(cnt[n_ter >= k, .(site, n_ter)], vt, by = "site"), file.path(out, st), sep = "\t")
  }
}
cat("done\n")
