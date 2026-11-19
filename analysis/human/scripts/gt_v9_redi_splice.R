#!/usr/bin/env Rscript
# gt_v9_redi_splice.R
# REDIportal v9 scoring: 58 published (_parts/pub) + 2 Cellsnp-lite hybrid (_parts/ch) configurations, as gt_v8_splice60.R
#   GT_V9_DIR/rediportal_v9/_parts/{pub,ch}/scoring_mismatch -> GT_V9_DIR/rediportal_v9/scoring_mismatch
# argument splice76 (2026-11-13): these 60 + 16 UMI-aware MarkDuplicates configurations (umimd_redi_score.sh)
#   GT_V9_DIR/rediportal_v9/{scoring_mismatch,umimd/scoring_mismatch} -> GT_V9_DIR/rediportal_v9/scoring_mismatch76
suppressPackageStartupMessages(library(data.table))
V9 = Sys.getenv("GT_V9_DIR", file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v9"))
R = file.path(V9, "rediportal_v9")
P76 = identical(commandArgs(TRUE)[1], "splice76")
a = file.path(R, "_parts/pub/scoring_mismatch"); b = file.path(R, "_parts/ch/scoring_mismatch"); o = file.path(R, "scoring_mismatch")
NB = 2L
if (P76) { a = o; b = file.path(R, "umimd/scoring_mismatch"); o = file.path(R, "scoring_mismatch76"); NB = 16L }
dir.create(o, recursive = TRUE, showWarnings = FALSE)
stopifnot(identical(unname(tools::md5sum(file.path(a, "positives.tsv"))), unname(tools::md5sum(file.path(b, "positives.tsv")))))
k = c("caller", "bam_type", "strand_mode", "snp_subtract")
for (f in list.files(a)) {
  if (grepl("^gt_tier_scoring_.*\\.tsv$|^unit_scores\\.tsv$", f)) {
    x = fread(file.path(a, f)); y = fread(file.path(b, f))
    stopifnot(identical(names(x), names(y)), nrow(merge(unique(x[, ..k]), unique(y[, ..k]))) == 0L, nrow(unique(y[, ..k])) == NB)
    fwrite(rbind(x, y), file.path(o, f), sep = "\t")
  } else file.copy(file.path(a, f), file.path(o, f), overwrite = TRUE)
}
if (P76) stopifnot(all(fread(file.path(b, "gt_tier_scoring_alt3.tsv"))$bam_type == "sampleUmiMd"))
cat("REDIportal v9 spliced:", nrow(unique(fread(file.path(o, "gt_tier_scoring_alt3.tsv"))[, ..k])), "configurations\n")
