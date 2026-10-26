#!/usr/bin/env Rscript
# =============================================================================
# gt_vartype_counts_v6.R
# Unique sites per mismatch type (varType, after stranding) in every NTC call file
# of the three local-mismatch-context states, all tiers:
#   none     = step 4 (cov/edit filtered)          gt_eval_v6/prewindow/results_<tier>/.../*_step4_cov_edit_filtered.tsv.gz
#   mismatch = step 5 (near other MM types removed) gt_eval_v6/prewindow/results_<tier>/.../*_step5_nomm.tsv.gz
#   both     = step 6 (+ clustered; benchmark set)  results_gt*_v6/.../*_pre100.tsv.gz
#   clustered = clustered only (pass 1 off, w2 searched)  results_gt*_v6_nomm/.../*_pre100.tsv.gz
#               (nomm_rerun_v6.sh, HANDOFF 2026-10-01 UPDATE-1)
# Output: gt_eval_v6/prewindow/vartype/vartype_counts_all_tiers.tsv.gz
#   window_filter, tier, unit, caller, bam_type, strand_mode, snp_subtract, varType, n
# Usage: Rscript gt_vartype_counts_v6.R [threads]
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
OUT <- file.path(PM, "gt_eval_v6", "prewindow", "vartype", "vartype_counts_all_tiers.tsv.gz")
argv <- commandArgs(trailingOnly = TRUE)
threads <- if (length(argv)) as.integer(argv[1]) else 1L

ROOT <- list(none     = c(arm = "gt_eval_v6/prewindow/results_arm", sample = "gt_eval_v6/prewindow/results_sample",
                          celltype = "gt_eval_v6/prewindow/results_celltype"),
             mismatch = c(arm = "gt_eval_v6/prewindow/results_arm", sample = "gt_eval_v6/prewindow/results_sample",
                          celltype = "gt_eval_v6/prewindow/results_celltype"),
             both     = c(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6"),
             clustered = c(arm = "results_gt_v6_nomm", sample = "results_gt_persample_v6_nomm",
                           celltype = "results_gt_celltype_v6_nomm"))
CALL_RE <- c(none = "_step4_cov_edit_filtered\\.tsv\\.gz$", mismatch = "_step5_nomm\\.tsv\\.gz$",
             both = "_pre100\\.tsv\\.gz$", clustered = "_pre100\\.tsv\\.gz$")
UNIT_RE <- c(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")

files <- rbindlist(lapply(names(ROOT), function(w) rbindlist(lapply(names(ROOT[[w]]), function(tr) {
  root <- file.path(PM, ROOT[[w]][[tr]])
  units <- grep(UNIT_RE[[tr]], list.dirs(root, recursive = FALSE, full.names = FALSE), value = TRUE)
  rbindlist(lapply(units, function(u) {
    base <- file.path(root, u, "04_filtered")
    f <- list.files(base, CALL_RE[[w]], recursive = TRUE, full.names = TRUE)
    if (!length(f)) return(NULL)
    p <- tstrsplit(sub(paste0("^", base, "/"), "", f), "/", fixed = TRUE)
    data.table(window_filter = w, tier = tr, unit = u, caller = p[[1]], bam_type = p[[2]],
               strand_mode = p[[3]], snp_subtract = sub("^snp_", "", p[[4]]), file = f)
  }))
}))))
message(nrow(files), " call files")
print(files[, .N, by = .(window_filter, tier)])

cnt <- rbindlist(mclapply(seq_len(nrow(files)), function(i) {
  f <- files$file[i]
  if (file.size(f) == 0) return(NULL)
  d <- unique(fread(f, select = c("site", "varType"), showProgress = FALSE))
  if (!nrow(d)) return(NULL)
  cbind(files[i, !"file"], d[, .(n = .N), by = varType])
}, mc.cores = threads))
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(cnt, OUT, sep = "\t")
message("wrote ", OUT, ": ", nrow(cnt), " rows")
