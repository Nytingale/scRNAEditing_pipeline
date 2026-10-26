#!/usr/bin/env Rscript
# ============================================================================
# WASP_refbias_source_data.R   (split from WASP_refbias_plot.R 2026-10-25)
#
# Source data of the per-dataset WASP reference-bias panels, at each coverage floor (10 / 20 / 50):
#   a_sites       per-site editing level with all reads and with WASP-passed reads, and dEL = EL_pass - EL_all,
#                 one row per site and floor (gzip)                          (WASP_hist_delta_editing_level)
#   a_stats       mean and median dEL and number of sites per floor
#   b_alleles     share of reads removed by WASP at REDIportal sites, by allele  (WASP_bar_removed_read_alleles)
#   b_gaps        reference minus edited removal gap (pp) with its 95% CI, per floor
#   c_sites       REDIportal sites covered without WASP vs still evaluable after it  (WASP_bar_site_attrition)
# WASP_panel_refbias stacks the three panels.
#
# INPUT : data/<ds>/WASP_persite_evaluable.tsv.gz, report/<ds>/WASP_{by_threshold,removed_alleles}.tsv
#         (WASP_prep_source_data.R)
# OUTPUT: source_data/<ds>/WASP_refbias.<part>.tsv (drawn by WASP_refbias_plot.R)
# Usage : Rscript WASP_refbias_source_data.R [dataset]   (default jzlxv3; datasets in WASP_datasets.tsv)
# ============================================================================
suppressPackageStartupMessages(library(data.table))
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(SCRIPT_DIR, "WASP_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
args <- commandArgs(trailingOnly = TRUE)
wasp_use_dataset(if (length(args) >= 1) args[1] else "jzlxv3")
STEM <- "WASP_refbias"

need <- function(f) {
  if (!file.exists(f)) stop("Missing ", f, "\nRun: Rscript ", file.path("scripts", "WASP_prep_source_data.R"))
  f
}
cov_facet <- function(th) {
  factor(sprintf("Coverage >= %d", th),
         levels = sprintf("Coverage >= %d", sort(unique(th))))
}

thr <- fread(need(file.path(REPORT_DIR, "WASP_by_threshold.tsv")))
THRESHOLDS <- sort(unique(thr$min_cov))

# (A) The evaluable table is written at the LOWEST floor and the floors nest, so each facet is a filter on it.
ev <- fread(need(file.path(DATA_DIR, "WASP_persite_evaluable.tsv.gz")))
dat <- rbindlist(lapply(THRESHOLDS, function(th) {
  ev[cov_all >= th & cov_pass >= th,
     .(min_cov = th, chr, pos, ref, EL_all, EL_pass, dEL)]
}))
dat[, cov_facet := cov_facet(min_cov)]
write_sd(dat, SRC_DIR, STEM, "a_sites", gzip = TRUE)
write_sd(dat[, .(avg = mean(dEL), med = median(dEL), n = .N), by = .(min_cov, cov_facet)], SRC_DIR, STEM, "a_stats")

# (B) Title case to match the cmap_allele keys.
al <- fread(need(file.path(REPORT_DIR, "WASP_removed_alleles.tsv")))
al[, allele := factor(stringr::str_to_title(allele), levels = c("Reference (A)", "Edited (G)"))]
al[, cov_facet := cov_facet(min_cov)]
write_sd(al, SRC_DIR, STEM, "b_alleles")
write_sd(thr[, .(min_cov, gap_pp, gap_ci_lo_pp, gap_ci_hi_pp)], SRC_DIR, STEM, "b_gaps")

# (C)
st <- melt(thr[, .(min_cov, Covered = n_sites_covered, Evaluable = n_sites_evaluable)],
           id.vars = "min_cov", variable.name = "set", value.name = "n_sites")
st[, set := factor(set, levels = c("Covered", "Evaluable"))]
st[, cov_facet := cov_facet(min_cov)]
st <- merge(st, thr[, .(min_cov, pct_sites_dropped)], by = "min_cov")
write_sd(st, SRC_DIR, STEM, "c_sites")
message("source data -> ", SRC_DIR)
