#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_ag_percent.R   (2026-10-13)
#
# A-to-G percentage per NTC unit and configuration on the v8 call state (mismatch = nearby-mismatch
# pass only: results_gt_{v8,persample_v8,celltype_v8}/<NTC unit>/04_filtered/**/*_step5_nomm.tsv.gz),
# for TERATOMA_v8_overall_vs_ag_plot.R. Same rule as ag_percent_units.tsv in
# gt_teratoma_plot_inputs_v6.R: unique site x varType per call file; an empty file counts as 0 sites.
# Output: gt_eval_v8/scoring/mismatch/ag_percent_units.tsv (tier, unit, config, n_sites, n_ag, ag_percent)
# Cell types: GT_UNIT_KEEP (default gt_eval_v8/unit_cov/celltype_keep_min20.txt; empty = all NTC cell types)
# Usage : Rscript gt_v8_ag_percent.R   (login node OK, ~2 min)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
PUB = c(arm = "results_gt_v8", sample = "results_gt_persample_v8", celltype = "results_gt_celltype_v8")
OUT = file.path(PM, "gt_eval_v8/scoring/mismatch/ag_percent_units.tsv")
KEEP = Sys.getenv("GT_UNIT_KEEP", file.path(PM, "gt_eval_v8/unit_cov/celltype_keep_min20.txt"))

ag = rbindlist(lapply(c("arm", "sample", "celltype"), function(tier) {
  root = file.path(PM, PUB[[tier]])
  units = grep("NTC_ctrl", list.dirs(root, recursive = FALSE, full.names = FALSE), value = TRUE)
  if (tier == "celltype" && nzchar(KEEP)) { keep = readLines(KEEP); stopifnot(all(keep %in% units)); units = keep }
  files = unlist(lapply(units, function(u)
    list.files(file.path(root, u, "04_filtered"), pattern = "_step5_nomm\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)))
  rbindlist(lapply(files, function(f) {
    # <unit>/04_filtered/<caller>/<bam_type>/<strand_mode>/snp_<on|off>/<file>
    p = strsplit(sub(paste0("^", root, "/"), "", f), "/", fixed = TRUE)[[1]]
    cfg = paste(p[3], p[4], p[5], sub("^snp_", "", p[6]), sep = "__")
    d = if (file.size(f) == 0) NULL else unique(fread(f, select = c("site", "varType"), showProgress = FALSE))
    data.table(tier = tier, unit = p[1], config = cfg, n_sites = NROW(d), n_ag = if (is.null(d)) 0L else sum(d$varType == "AG"))
  }))
}))
stopifnot(ag[, .N, by = .(tier, unit, config)][, all(N == 1L)])
ag[, ag_percent := fifelse(n_sites > 0, 100 * n_ag / n_sites, NA_real_)]
fwrite(ag, OUT, sep = "\t")
cat(nrow(ag), "call files;", uniqueN(ag$config), "configurations; median A-to-G % by tier:",
    ag[, .(m = round(median(ag_percent, na.rm = TRUE), 1)), by = tier][, paste(tier, m)], "\n")
