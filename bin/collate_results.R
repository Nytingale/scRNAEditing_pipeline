#!/usr/bin/env Rscript
#
# collate_results.R — RESULTS GENERATION ONLY. NO PLOTTING HAPPENS HERE.
#
# Writes the benchmarking summary table and nothing else. Figures are built
# separately, by the publication deliverable:
#     Rscript publication/benchmark_grid/scripts/BENCHMARK_prep_source_data.R
#     bash    publication/benchmark_grid/scripts/run_BENCHMARK_figures.sh
#
# WHY THEY ARE SEPARATE. This script used to system2() the plotting script, so
# plotting ran inside COLLATE_RESULTS, the pipeline's final step. A plotting bug
# therefore failed the whole pipeline: a 35-second ggplot OOM discarded ~28 h of
# compute and cost two full relaunches of nf_pm_ps (jobs 1693635, 1705676) even
# though every result had already been computed. Figures must be regenerable
# without re-running the pipeline, and a broken figure must never be able to
# fail a multi-day run.
#
options(scipen = 999)
suppressPackageStartupMessages({
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: collate_results.R <all_run_metrics.tsv> <output_dir>")

metrics_file = args[1]
outdir       = args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dt = fread(metrics_file)
if (nrow(dt) == 0) {
  message("No runs in metrics file. Nothing to collate.")
  fwrite(data.table(), file.path(outdir, "all_results_summary.tsv"), sep = "\t")
  quit(status = 0)
}

dt = dt[sites_step1_multiMM > 0]

# `snp_subtract` is part of the RUN KEY, not an incidental parameter: the grid
# runs each (sample, caller, bam_type, strand_mode) both with and without SNP
# subtraction. Dropping it from the summary left two rows per key that were
# identical in every retained column, which made
# (sample, caller, bam_type[, strand_mode]) non-unique. Any downstream merge on
# those columns then produced a CARTESIAN product -- the strand-mode comparison
# panels turned 220 real comparisons into 880 points, half of them cross-pairing
# annotation-with-subtraction against hybrid-without (not a meaningful
# comparison), and the inflated count also pushed the panel past
# raster_point_min into a 292 M px raster that OOM-killed COLLATE_RESULTS twice.
key_cols = c("sample", "caller", "bam_type", "strand_mode", "snp_subtract",
             "best_w1", "best_w2",
             "sites_step1_multiMM", "sites_step2_noSNP", "sites_step3_stranded",
             "sites_step4_covEdit", "sites_step5_noMM", "sites_step6_clustered",
             "sites_step7_no100pct",
             "ag_pct_step1", "ag_pct_step2", "ag_pct_step3", "ag_pct_step4",
             "ag_pct_step5", "ag_pct_step6", "ag_pct_step7",
             "ag_n_step7",
             "n_strand_from_read", "n_strand_from_annot", "n_strand_from_caller")
present_cols = intersect(key_cols, colnames(dt))
summary_file = file.path(outdir, "all_results_summary.tsv")
fwrite(dt[, ..present_cols], summary_file, sep = "\t", quote = FALSE)

# The run key must stay unique in the summary, because every downstream consumer
# joins or groups on it. A duplicate here is what produced the Cartesian merge
# described above, so it is checked at the point the table is written rather than
# discovered later in a figure.
grid_cols = intersect(c("sample", "caller", "bam_type", "strand_mode", "snp_subtract",
                        "min_cov", "min_alt", "min_edit_level",
                        "jacusa2_min_score", "remove_100pct"), present_cols)
if (length(grid_cols)) {
  n_key = nrow(unique(dt[, ..grid_cols]))
  if (n_key < nrow(dt)) {
    warning("run key is NOT unique in the summary: ", nrow(dt), " rows collapse to ",
            n_key, " keys on (", paste(grid_cols, collapse = ", "),
            "). A grid dimension is missing from key_cols; any downstream merge ",
            "on it will produce a Cartesian product.")
  }
  message("run key  : ", paste(grid_cols, collapse = ", "))
  message("unique   : ", n_key, " of ", nrow(dt), " rows")
}

message("Wrote summary: ", summary_file, " (", nrow(dt), " runs)")
message("Collation complete. NO figures were produced here — build them with:")
message("  Rscript publication/benchmark_grid/scripts/BENCHMARK_prep_source_data.R ", summary_file)
message("  bash    publication/benchmark_grid/scripts/run_BENCHMARK_figures.sh")
