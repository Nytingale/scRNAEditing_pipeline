#!/usr/bin/env Rscript
# =============================================================================
# gt_vartype_counts_v8.R   (2026-10-15; from gt_vartype_counts_v6.R)
# Unique sites per mismatch type (varType, after stranding) in every NTC call file of the four
# local-mismatch-context states, v8 caller runs (2026-10-14; BCFtools BAQ on, REF=N fixed):
#   none      = step 4      results_gt{,_persample,_celltype}_v8/.../*_step4_cov_edit_filtered.tsv.gz
#   mismatch  = step 5      results_gt{,_persample,_celltype}_v8/.../*_step5_nomm.tsv.gz
#   both      = step 6      results_gt{,_persample,_celltype}_v8/.../*_pre100.tsv.gz
#   clustered = pass 1 off  results_gt{,_persample,_celltype}_v8_nomm/.../*_pre100.tsv.gz (nomm_rerun_v8.sh)
# Cell types: the 11 of gt_eval_v8/unit_cov/celltype_keep_min20.txt (as the v8 scoring)
# Output: gt_eval_v8/vartype/vartype_counts_all_tiers.tsv.gz
#   window_filter, tier, unit, caller, bam_type, strand_mode, snp_subtract, varType, n
# Usage: Rscript gt_vartype_counts_v8.R [threads]
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
OUT <- file.path(PM, "gt_eval_v8", "vartype", "vartype_counts_all_tiers.tsv.gz")
KEEP <- readLines(file.path(PM, "gt_eval_v8", "unit_cov", "celltype_keep_min20.txt"))
argv <- commandArgs(trailingOnly = TRUE)
threads <- if (length(argv)) as.integer(argv[1]) else 1L

PUB <- c(arm = "results_gt_v8", sample = "results_gt_persample_v8", celltype = "results_gt_celltype_v8")
ROOT <- list(none = PUB, mismatch = PUB, both = PUB, clustered = setNames(paste0(PUB, "_nomm"), names(PUB)))
CALL_RE <- c(none = "_step4_cov_edit_filtered\\.tsv\\.gz$", mismatch = "_step5_nomm\\.tsv\\.gz$",
             both = "_pre100\\.tsv\\.gz$", clustered = "_pre100\\.tsv\\.gz$")
UNIT_RE <- c(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")

files <- rbindlist(lapply(names(ROOT), function(w) rbindlist(lapply(names(ROOT[[w]]), function(tr) {
  root <- file.path(PM, ROOT[[w]][[tr]])
  units <- grep(UNIT_RE[[tr]], list.dirs(root, recursive = FALSE, full.names = FALSE), value = TRUE)
  if (tr == "celltype") { stopifnot(all(KEEP %in% units)); units <- KEEP }
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
print(dcast(files[, .N, by = .(window_filter, tier)], tier ~ window_filter, value.var = "N"))
stopifnot(files[, .N, by = .(window_filter, tier)][, uniqueN(N), by = tier][, all(V1 == 1L)])

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
