#!/usr/bin/env Rscript
# ============================================================================
# BENCHMARK_prep_source_data.R — source-data generation for the benchmark-grid
# deliverable. NO PLOTTING HAPPENS HERE.
# ============================================================================
# Stage 2 of 3:
#   Stage 1  the Nextflow pipeline        FASTQ/BAM -> per-run filtering metrics
#                                         (summary/all_run_metrics.tsv, then
#                                         COLLATE_RESULTS -> all_results_summary.tsv)
#   Stage 2  THIS SCRIPT                  per-run table -> sample-aggregated
#                                         plot-ready tables + report tables
#   Stage 3  run_BENCHMARK_figures.sh     tables -> panels
#
# WHAT THIS DELIVERABLE COMPARES
# The pipeline runs a full grid of variant caller x BAM preprocessing x strand
# assignment x SNP subtraction, over every sample, and records site counts and
# the A-to-G fraction at each of 7 filtering steps. A-to-G fraction is the
# specificity readout: real A-to-I editing appears as A->G, so the closer to
# 100% A-to-G, the less the surviving call set is contaminated by other
# mismatch types.
#
# >>> AGGREGATION OVER `sample` IS THE DESIGN, AND IS RECORDED. <<<
# `sample` is the biological REPLICATE dimension, so every output here is a
# median with an interquartile range across samples, carrying `n_samples`. Every
# OTHER grid dimension (caller, bam_type, strand_mode, snp_subtract) is kept as a
# separate row, because collapsing one of those would compare different methods
# as if they were one.
# The predecessor of this deliverable plotted all 1280 runs at once, which needed
# 120 facet columns in 18 inches (0.15 in per facet) and was unreadable; it also
# dropped `snp_subtract` entirely, which fabricated line trajectories, hid bars
# behind bars, and summed stacked bars to double height.
#
# Usage:
#   Rscript BENCHMARK_prep_source_data.R [DATASET] [SUMMARY_TSV]
# DATASET defaults to pm_gt_persample and must exist in
# scripts/BENCHMARK_datasets.tsv, which also supplies its summary table path.
# SUMMARY_TSV overrides the registry for a one-off run.
#
# Reads:
#   <SUMMARY_TSV>                        one row per run, written by COLLATE_RESULTS
#
# Writes (data/, plot-ready):
#   BENCHMARK_final_by_run.tsv.gz        per-run final-step metrics (the rows the
#                                        box panel plots; sample kept, since a
#                                        box IS the across-sample distribution)
#   BENCHMARK_attrition_by_step.tsv.gz   median/IQR of site count and A-to-G %
#                                        per grid cell x filtering step
#   BENCHMARK_strand_mode_cmp.tsv        annotation vs hybrid, paired per run
#                                        then aggregated over samples
#   BENCHMARK_windows.tsv                adaptive window (w1,w2) per grid cell
#   BENCHMARK_strand_source.tsv          hybrid strand-evidence source split
# Writes (report/):
#   BENCHMARK_summary.tsv                headline metrics, long format
#   BENCHMARK_grid_scope.tsv             the grid actually present in the input
#   BENCHMARK_data_dictionary.tsv        column definitions for all of the above
# Writes (source_data/):
#   SOURCE_DATA_MANIFEST.tsv             primary inputs + which stage made them
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

script_dir <- (function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
})()
DELIV_DIR <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
REPO_DIR  <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

# Dataset and the four per-dataset output trees come from the registry via
# BENCHMARK_common.R, so this stage and the plotting stage cannot disagree about
# where a dataset's tables live. Previously these were FLAT, so a second dataset
# silently overwrote the first.
source(file.path(script_dir, "BENCHMARK_common.R"))
DATASET <- if (length(args) >= 1) args[1] else "pm_gt_persample"
bench_set_dataset(DATASET)
# The registry supplies the input table; an explicit second argument overrides it
# (useful for a one-off run against a non-registered summary).
if (length(args) >= 2) SUMMARY_TSV <- args[2]

if (!file.exists(SUMMARY_TSV)) {
  stop("summary table not found: ", SUMMARY_TSV,
       "\nRun the pipeline (COLLATE_RESULTS writes it), or fix the path in ",
       basename(BENCHMARK_DATASETS))
}
message("dataset : ", DATASET, "  (", DATASET_LABEL, ")")
message("input   : ", SUMMARY_TSV)

dt <- fread(SUMMARY_TSV)
if ("sites_step1_multiMM" %in% names(dt)) dt <- dt[sites_step1_multiMM > 0]
if (!nrow(dt)) stop("no non-empty runs in ", SUMMARY_TSV)

GRID_DIMS <- c("caller", "bam_type", "strand_mode", "snp_subtract")
AGG_DIMS  <- "sample"
GRID      <- intersect(GRID_DIMS, names(dt))
STEP_LABELS <- c("1:multiMM", "2:noSNP", "3:stranded", "4:covEdit",
                 "5:noMM(w1)", "6:clustered(w2)", "7:no100pct")

# ── Grid scope: what is actually in this input ───────────────────────────────
# >>> The completeness check the project rule demands. <<<
# Recorded BEFORE anything is aggregated, so the report states exactly which
# combinations exist, and a missing dimension is caught here rather than in a
# figure. A non-unique run key means a downstream merge would go Cartesian.
missing_dims <- setdiff(GRID_DIMS, names(dt))
if (length(missing_dims)) {
  warning("grid dimension(s) absent from the summary: ", paste(missing_dims, collapse = ", "),
          " — combinations differing only in those are indistinguishable here")
}
key_cols <- c(AGG_DIMS, GRID)
n_key <- nrow(unique(dt[, ..key_cols]))
if (n_key < nrow(dt)) {
  stop("run key (", paste(key_cols, collapse = ", "), ") is NOT unique: ",
       nrow(dt), " rows collapse to ", n_key,
       " keys. Aggregating would silently average unrelated runs.")
}
message("run key : ", paste(key_cols, collapse = ", "), "  (", n_key, " unique of ",
        nrow(dt), " rows)")

scope <- rbindlist(lapply(c(AGG_DIMS, GRID), function(cl) data.table(
  dataset = DATASET, dimension = cl, role = if (cl %in% AGG_DIMS) "aggregated" else "distinguished",
  n_levels = uniqueN(dt[[cl]]),
  levels = paste(sort(unique(as.character(dt[[cl]]))), collapse = "|"))))
scope <- rbindlist(list(scope, data.table(
  dataset = DATASET, dimension = "(runs)", role = "total",
  n_levels = nrow(dt), levels = NA_character_)), use.names = TRUE)
fwrite(scope, file.path(REPORT_DIR, "BENCHMARK_grid_scope.tsv"), sep = "\t")

# ── Helper: median / IQR over the replicate dimension ───────────────────────
agg_over_samples <- function(d, value, by) {
  d <- d[!is.na(get(value))]
  d[, .(n_samples = uniqueN(sample),
        median = as.numeric(median(get(value))),
        q25    = as.numeric(quantile(get(value), 0.25, names = FALSE)),
        q75    = as.numeric(quantile(get(value), 0.75, names = FALSE)),
        mean   = mean(get(value)),
        sd     = if (.N > 1L) sd(get(value)) else NA_real_),
    by = by]
}

# ── data/: per-run final-step metrics (box panel input) ─────────────────────
# `sample` is KEPT here: a boxplot's whole purpose is to show the across-sample
# distribution, so the rows passed to ggplot are the per-sample values and the
# aggregation is done by the boxplot stat. Source data must be the exact rows
# plotted (publication/README.md), which is why this is not pre-aggregated.
final_cols <- intersect(c("sample", GRID, "ag_pct_step7", "sites_step7_no100pct",
                          "ag_n_step7", "best_w1", "best_w2"), names(dt))
final_by_run <- dt[, ..final_cols]
final_by_run[, ag_pct_step7_pct := 100 * ag_pct_step7]
fwrite(final_by_run, file.path(DATA_DIR, "BENCHMARK_final_by_run.tsv.gz"), sep = "\t")
message("wrote   : data/BENCHMARK_final_by_run.tsv.gz (", nrow(final_by_run), " runs)")

# ── data/: attrition + A-to-G progression, median/IQR per grid cell x step ──
site_cols <- grep("^sites_step", names(dt), value = TRUE)
ag_cols   <- grep("^ag_pct_step", names(dt), value = TRUE)

long_sites <- melt(dt, id.vars = c("sample", GRID), measure.vars = site_cols,
                   variable.name = "step_col", value.name = "n_sites")
long_sites[, step_num := as.integer(sub("sites_step([0-9]+)_.*", "\\1", step_col))]

long_ag <- melt(dt, id.vars = c("sample", GRID), measure.vars = ag_cols,
                variable.name = "step_col", value.name = "ag_pct")
long_ag[, step_num := as.integer(sub("ag_pct_step([0-9]+)", "\\1", step_col))]

a_sites <- agg_over_samples(long_sites, "n_sites", c(GRID, "step_num"))
set(a_sites, j = "metric", value = "n_sites")
a_ag <- agg_over_samples(long_ag, "ag_pct", c(GRID, "step_num"))
# set() rather than [[<- : the latter copies the whole table and then triggers
# data.table's shallow-copy warning on the next :=.
for (cl in c("median", "q25", "q75", "mean", "sd")) set(a_ag, j = cl, value = 100 * a_ag[[cl]])
set(a_ag, j = "metric", value = "ag_pct")

attrition <- rbindlist(list(a_sites, a_ag), use.names = TRUE)
attrition[, step_label := factor(STEP_LABELS[step_num], levels = STEP_LABELS)]
setcolorder(attrition, c(GRID, "metric", "step_num", "step_label",
                         "n_samples", "median", "q25", "q75", "mean", "sd"))
fwrite(attrition, file.path(DATA_DIR, "BENCHMARK_attrition_by_step.tsv.gz"), sep = "\t")
message("wrote   : data/BENCHMARK_attrition_by_step.tsv.gz (", nrow(attrition), " rows)")

# ── data/: annotation vs hybrid, PAIRED WITHIN a run, then aggregated ───────
# Pairing happens per (sample, caller, bam_type, snp_subtract) BEFORE any
# aggregation, so each difference compares two strand modes of the SAME run.
# snp_subtract is part of the pairing key: without it the merge is a Cartesian
# product that pairs annotation-with-subtraction against hybrid-without.
strand_cmp <- NULL
if ("strand_mode" %in% names(dt) && uniqueN(dt$strand_mode) > 1L) {
  pair_by <- intersect(c("sample", "caller", "bam_type", "snp_subtract"), names(dt))
  a <- dt[strand_mode == "annotation",
          c(mget(pair_by), .(ag_annot = 100 * ag_pct_step7, sites_annot = sites_step7_no100pct))]
  h <- dt[strand_mode == "hybrid",
          c(mget(pair_by), .(ag_hybrid = 100 * ag_pct_step7, sites_hybrid = sites_step7_no100pct))]
  if (nrow(a) && nrow(h)) {
    if (anyDuplicated(a[, ..pair_by]) || anyDuplicated(h[, ..pair_by]))
      stop("strand-mode pairing key is not unique; the merge would be Cartesian")
    pr <- merge(a, h, by = pair_by)
    pr[, `:=`(d_ag = ag_hybrid - ag_annot, d_sites = sites_hybrid - sites_annot)]
    by2 <- setdiff(pair_by, "sample")
    strand_cmp <- merge(
      agg_over_samples(pr, "d_ag", by2)[, .(n_samples, d_ag_median = median,
                                            d_ag_q25 = q25, d_ag_q75 = q75), by = by2],
      agg_over_samples(pr, "d_sites", by2)[, .(d_sites_median = median,
                                               d_sites_q25 = q25, d_sites_q75 = q75), by = by2],
      by = by2)
    fwrite(pr, file.path(DATA_DIR, "BENCHMARK_strand_mode_paired_by_run.tsv.gz"), sep = "\t")
    fwrite(strand_cmp, file.path(DATA_DIR, "BENCHMARK_strand_mode_cmp.tsv"), sep = "\t")
    message("wrote   : data/BENCHMARK_strand_mode_cmp.tsv (", nrow(strand_cmp),
            " grid cells from ", nrow(pr), " paired runs)")
  }
}

# ── data/: adaptive windows ─────────────────────────────────────────────────
# w1/w2 are fitted per run, so the across-sample summary is the MEDIAN window.
if (all(c("best_w1", "best_w2") %in% names(dt))) {
  win <- dt[, .(n_samples = uniqueN(sample),
                w1_median = as.numeric(median(best_w1)),
                w2_median = as.numeric(median(best_w2)),
                w1_min = min(best_w1), w1_max = max(best_w1),
                w2_min = min(best_w2), w2_max = max(best_w2),
                ag_pct_step6_median = 100 * as.numeric(median(ag_pct_step6))), by = GRID]
  fwrite(win, file.path(DATA_DIR, "BENCHMARK_windows.tsv"), sep = "\t")
  message("wrote   : data/BENCHMARK_windows.tsv (", nrow(win), " grid cells)")
}

# ── data/: hybrid strand-evidence source ────────────────────────────────────
if (all(c("n_strand_from_read", "n_strand_from_annot") %in% names(dt))) {
  hy <- dt[strand_mode == "hybrid" & !is.na(n_strand_from_read)]
  if (nrow(hy)) {
    ss <- melt(hy, id.vars = c("sample", GRID),
               measure.vars = c("n_strand_from_read", "n_strand_from_annot"),
               variable.name = "source", value.name = "n_sites")
    ss[, source := sub("n_strand_from_", "", source)]
    ss_agg <- agg_over_samples(ss, "n_sites", c(setdiff(GRID, "strand_mode"), "source"))
    fwrite(ss_agg, file.path(DATA_DIR, "BENCHMARK_strand_source.tsv"), sep = "\t")
    message("wrote   : data/BENCHMARK_strand_source.tsv (", nrow(ss_agg), " rows)")
  }
}

# ── report/: headline metrics ───────────────────────────────────────────────
best <- final_by_run[, .(ag_pct_median = median(100 * ag_pct_step7),
                         sites_median = median(sites_step7_no100pct),
                         n_samples = uniqueN(sample)), by = GRID][order(-ag_pct_median)]
fwrite(best, file.path(REPORT_DIR, "BENCHMARK_final_by_grid_cell.tsv"), sep = "\t")

sm <- function(metric, value, note = NA_character_)
  data.table(dataset = DATASET, metric = metric, value = as.character(value), note = note)
summary_rows <- rbindlist(list(
  sm("n_runs", nrow(dt), "rows in the input summary"),
  sm("n_samples", uniqueN(dt$sample), "replicate dimension, aggregated over"),
  sm("n_grid_cells", nrow(unique(dt[, ..GRID])), paste(GRID, collapse = " x ")),
  sm("n_callers", uniqueN(dt$caller)),
  sm("n_bam_types", uniqueN(dt$bam_type)),
  sm("best_cell_ag_pct", sprintf("%.2f", best$ag_pct_median[1]),
     paste0("highest median final A-to-G %: ",
            paste(unlist(best[1, ..GRID]), collapse = " / "))),
  sm("worst_cell_ag_pct", sprintf("%.2f", best$ag_pct_median[nrow(best)]),
     paste0("lowest median final A-to-G %: ",
            paste(unlist(best[nrow(best), ..GRID]), collapse = " / ")))
), use.names = TRUE)
if (!is.null(strand_cmp)) {
  summary_rows <- rbindlist(list(summary_rows,
    sm("strand_hybrid_minus_annot_ag_median",
       sprintf("%+.3f", median(strand_cmp$d_ag_median)),
       "median across grid cells of the paired hybrid - annotation A-to-G % shift")),
    use.names = TRUE)
}
fwrite(summary_rows, file.path(REPORT_DIR, "BENCHMARK_summary.tsv"), sep = "\t")
message("wrote   : report/BENCHMARK_summary.tsv")

# ── report/: data dictionary ────────────────────────────────────────────────
dict <- rbindlist(list(
  data.table(file = "all (grid columns)", column = "caller",       meaning = "Variant caller: bcftools, gatk, jacusa2, reditools3, scallele."),
  data.table(file = "all (grid columns)", column = "bam_type",     meaning = "BAM preprocessing variant (sample, sampleRmDup, sampleUmiDup, recal, recalq20, recaldefault)."),
  data.table(file = "all (grid columns)", column = "strand_mode",  meaning = "Strand assignment: annotation-only, or hybrid (read-strand evidence with annotation fallback)."),
  data.table(file = "all (grid columns)", column = "snp_subtract", meaning = "Whether known-SNP subtraction was applied (on/off). A GRID DIMENSION, not a fixed setting."),
  data.table(file = "all (aggregated)",   column = "n_samples",    meaning = "Number of samples contributing to the median/IQR. The replicate count."),
  data.table(file = "all (aggregated)",   column = "median,q25,q75,mean,sd", meaning = "Across-SAMPLE summary of the metric for that grid cell. Samples are the replicate dimension."),
  data.table(file = "BENCHMARK_attrition_by_step.tsv.gz", column = "metric",     meaning = "'n_sites' = surviving site count; 'ag_pct' = A-to-G percentage of mismatches."),
  data.table(file = "BENCHMARK_attrition_by_step.tsv.gz", column = "step_num,step_label", meaning = "Filtering step 1-7 and its label; the pipeline applies them in this order."),
  data.table(file = "BENCHMARK_final_by_run.tsv.gz", column = "ag_pct_step7,ag_pct_step7_pct", meaning = "Final A-to-G fraction and the same value as a percentage. The specificity readout."),
  data.table(file = "BENCHMARK_final_by_run.tsv.gz", column = "sites_step7_no100pct", meaning = "Sites surviving all 7 filters, after removing 100%-edited sites."),
  data.table(file = "BENCHMARK_final_by_run.tsv.gz", column = "sample",      meaning = "KEPT here on purpose: the box panel plots per-sample values and lets the boxplot stat do the aggregating."),
  data.table(file = "BENCHMARK_strand_mode_paired_by_run.tsv.gz", column = "d_ag,d_sites", meaning = "hybrid minus annotation, PAIRED WITHIN one run. Positive d_ag = hybrid raises the A-to-G fraction."),
  data.table(file = "BENCHMARK_windows.tsv", column = "w1_median,w2_median", meaning = "Median fitted adaptive-window parameters (w1 = step 5 mismatch window, w2 = step 6 cluster window)."),
  data.table(file = "BENCHMARK_strand_source.tsv", column = "source",   meaning = "Origin of the strand call in hybrid mode: 'read' (read-strand evidence) or 'annot' (annotation fallback)."),
  data.table(file = "BENCHMARK_grid_scope.tsv", column = "role",        meaning = "'distinguished' = kept as its own row/facet; 'aggregated' = deliberately collapsed (the replicate dimension).")
), use.names = TRUE)
dict[, dataset := DATASET]
fwrite(dict, file.path(REPORT_DIR, "BENCHMARK_data_dictionary.tsv"), sep = "\t")
message("wrote   : report/BENCHMARK_data_dictionary.tsv")

# ── source_data/: manifest of primary inputs ────────────────────────────────
man <- data.table(
  dataset = DATASET,
  file = c(SUMMARY_TSV),
  produced_by = c("Nextflow COLLATE_RESULTS (nextflow/bin/collate_results.R)"),
  stage = c("1"),
  rows = c(nrow(dt)),
  note = c("one row per run; the run key is (sample, caller, bam_type, strand_mode, snp_subtract)")
)
fwrite(man, file.path(SRC_DIR, "SOURCE_DATA_MANIFEST.tsv"), sep = "\t")
message("wrote   : source_data/SOURCE_DATA_MANIFEST.tsv")

message("\nStage 2 complete. Build the figures with:")
message("  bash ", file.path("publication/benchmark_grid/scripts", "run_BENCHMARK_figures.sh"))
