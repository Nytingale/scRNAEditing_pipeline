#!/usr/bin/env Rscript
# gt_v9_filter_sweep_ag_pct.R   (2026-11-16)
# Call-set size and A-to-G % of each of the 76 configurations in the merged NTC arm under each downstream-filter
# setting of the filter sweep (filter_sweep76_v9.sh), from the unique (site, varType) of each step5_nomm call set as in
# gt_v10_configs76.R; and, against the main filters (centre), the calls a setting adds and their A-to-G %.
# Input : human/pm/gt_eval_v9/filter_sweep76/results/<setting>/results_gt_v8/NTC_ctrl/04_filtered/
# Output: human/pm/gt_eval_v9/filter_sweep76/ag_pct_arm.tsv (setting, config, n_calls, ag_pct, n_added, ag_pct_added)
#         printed: per setting, median and range over the 76
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("SCRNA_DIR", "."))
W <- "human/pm/gt_eval_v9/filter_sweep76"
SET <- c("centre", "cov5", "cov15", "alt2", "alt4", "edit2", "edit3", "window_off", "all_off")
rd <- function(s) {
  d <- file.path(W, "results", s, "results_gt_v8", "NTC_ctrl", "04_filtered")
  f <- list.files(d, "_step5_nomm\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  p <- tstrsplit(sub(paste0("^", d, "/"), "", dirname(f)), "/", fixed = TRUE)
  x <- data.table(file = f, config = paste(p[[1]], p[[2]], p[[3]], sub("^snp_", "", p[[4]]), sep = "__"))
  stopifnot(nrow(x) == 76L, !anyDuplicated(x$config))
  x[, {
    if (file.size(file) == 0) data.table(site = character(), varType = character())
    else unique(fread(file, select = c("site", "varType")))
  }, by = config][, setting := s]
}
SET <- SET[dir.exists(file.path(W, "results", SET))]
calls <- rbindlist(lapply(SET, rd))
cen <- calls[setting == "centre", .(config, site, varType, in_centre = TRUE)]
calls <- merge(calls, cen, by = c("config", "site", "varType"), all.x = TRUE)
calls[is.na(in_centre), in_centre := FALSE]
res <- calls[, .(n_calls = .N, ag_pct = 100 * mean(varType == "AG"),
                 n_added = sum(!in_centre), ag_pct_added = 100 * mean(varType[!in_centre] == "AG")), by = .(setting, config)]
fwrite(res, file.path(W, "ag_pct_arm.tsv"), sep = "\t")
options(width = 200)
print(res[, .(median_n = median(n_calls), median_ag_pct = round(median(ag_pct), 1),
              ag_pct_range = sprintf("%.1f-%.1f", min(ag_pct), max(ag_pct)),
              median_added = median(n_added), median_ag_pct_added = round(median(ag_pct_added, na.rm = TRUE), 1),
              ag_pct_added_range = sprintf("%.1f-%.1f", min(ag_pct_added, na.rm = TRUE), max(ag_pct_added, na.rm = TRUE))),
          by = setting])
d <- dcast(res[setting %in% c("centre", "all_off")], config ~ setting, value.var = c("ag_pct", "n_calls"))
cat(sprintf("all_off vs centre: A-to-G %% lower in %d of 76; call sets %.1f-%.1f-fold larger (median %.1f)\n",
            d[, sum(ag_pct_all_off < ag_pct_centre)], d[, min(n_calls_all_off / n_calls_centre)],
            d[, max(n_calls_all_off / n_calls_centre)], d[, median(n_calls_all_off / n_calls_centre)]))
