#!/usr/bin/env Rscript
# ============================================================================
# WASP_refbias_platforms_source_data.R   (split from WASP_refbias_plot_platforms.R 2026-10-25)
# Source data for response figure 15 (WASP reference-bias correction), both platforms and all coverage floors:
#   a_delta    per-site editing level shift, EL(WASP-passed) - EL(all reads), per platform x coverage floor
#   a_stats    mean, median and number of sites per facet
#   b_removed  reads removed at REDIportal sites, per allele, platform and coverage floor
#   b_gaps     reference - edited removal gap with its 95% CI clustered by site (WASP_gap_ci.sh)
#   c_sites    A-to-I sites covered vs evaluable after WASP
#   c_drops    sites lost per platform and coverage floor
# INPUT : data/<ds>/WASP_persite_evaluable.tsv.gz, report/<ds>/WASP_{summary,by_threshold,removed_alleles,gap_ci}.tsv
#         (WASP_prep_source_data.R, WASP_gap_ci.sh; ds = jzlxv3, jzlss2)
# OUTPUT: source_data/combined/WASP_refbias_platforms.<part>.tsv (a_delta gzipped)
#         (drawn by WASP_refbias_platforms_plot.R)
# Usage : Rscript WASP_refbias_platforms_source_data.R [dataset ...]   (default jzlxv3 jzlss2)
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "WASP_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
STEM <- "WASP_refbias_platforms"

args <- commandArgs(trailingOnly = TRUE)
DATASETS <- if (length(args)) args else c("jzlxv3", "jzlss2")
cfgs <- lapply(DATASETS, wasp_dataset)
names(cfgs) <- DATASETS
need <- function(f, what) {
  if (!file.exists(f)) stop("Missing ", f, "\nRun: ", what)
  f
}
plat_label <- function(cfg) sprintf("%s (%s)", sub("^Human PBMC, ", "", cfg$label), cfg$dataset)
PLAT_LEVELS <- vapply(cfgs, plat_label, character(1))
plat_facet <- function(ds) factor(PLAT_LEVELS[ds], levels = PLAT_LEVELS)
plat_short <- function(ds) sub(" \\(.*", "", PLAT_LEVELS[ds])

MINCOVS <- vapply(DATASETS, function(ds) {
  s = fread(need(file.path(cfgs[[ds]]$report_dir, "WASP_summary.tsv"),
                 paste("Rscript scripts/WASP_prep_source_data.R", ds)))
  as.integer(s[metric == "min_coverage", value])
}, integer(1))
if (length(unique(MINCOVS)) != 1L) {
  stop("Datasets were built at different coverage floors (",
       paste(sprintf("%s=%d", DATASETS, MINCOVS), collapse = ", "),
       "). A side-by-side figure would be comparing unlike sets.")
}
thr_all <- rbindlist(lapply(DATASETS, function(ds) {
  t = fread(need(file.path(cfgs[[ds]]$report_dir, "WASP_by_threshold.tsv"),
                 paste("Rscript scripts/WASP_prep_source_data.R", ds)))
  t[, dataset := ds][]
}))
THRESHOLDS <- sort(unique(thr_all$min_cov))
cov_facet <- function(th) factor(sprintf("Coverage >= %d", th), levels = sprintf("Coverage >= %d", THRESHOLDS))
wasp_use_combined()

# a: per-site editing level shift
dat.a <- rbindlist(lapply(DATASETS, function(ds) {
  ev = fread(need(file.path(cfgs[[ds]]$tab_dir, "WASP_persite_evaluable.tsv.gz"),
                  paste("Rscript scripts/WASP_prep_source_data.R", ds)))
  rbindlist(lapply(THRESHOLDS, function(th) {
    ev[cov_all >= th & cov_pass >= th,
       .(dataset = ds, min_cov = th, chr, pos, ref, EL_all, EL_pass, dEL)]
  }))
}))
dat.a[, `:=`(platform = plat_facet(dataset), cov = cov_facet(min_cov))]
stats.a <- dat.a[, .(avg = mean(dEL), med = median(dEL), n = .N), by = .(dataset, platform, min_cov, cov)]
stats.a[, lab := sprintf("Mean = %+.4f\nMedian = %+.4f\n%s sites", avg, med, format(n, big.mark = ","))]
write_sd(dat.a, SRC_DIR, STEM, "a_delta", gzip = TRUE)
write_sd(stats.a, SRC_DIR, STEM, "a_stats")

# b: read loss per allele; the gap is site-clustered (reads at one site are not independent)
dat.b <- rbindlist(lapply(DATASETS, function(ds) {
  d = fread(need(file.path(cfgs[[ds]]$report_dir, "WASP_removed_alleles.tsv"),
                 paste("Rscript scripts/WASP_prep_source_data.R", ds)))
  d[, dataset := ds][]
}))
dat.b[, allele := factor(stringr::str_to_title(allele), levels = c("Reference (A)", "Edited (G)"))]
dat.b[, `:=`(platform = plat_facet(dataset), cov = cov_facet(min_cov))]
gaps.b <- rbindlist(lapply(DATASETS, function(ds) {
  g = fread(need(file.path(cfgs[[ds]]$report_dir, "WASP_gap_ci.tsv"), paste("qsub scripts/WASP_gap_ci.sh", ds)))
  g[, dataset := ds][]
}))
gaps.b[, lab := sprintf("%s: %.2f pp [%.2f, %.2f]", plat_short(dataset), gap_pp, cluster_lo_pp, cluster_hi_pp)]
gaps.b <- gaps.b[, .(lab = paste(lab, collapse = "\n")), by = min_cov]
gaps.b[, cov := cov_facet(min_cov)]
write_sd(dat.b, SRC_DIR, STEM, "b_removed")
write_sd(gaps.b, SRC_DIR, STEM, "b_gaps")

# c: site attrition
dat.c <- melt(thr_all[, .(dataset, min_cov, Covered = n_sites_covered, Evaluable = n_sites_evaluable)],
              id.vars = c("dataset", "min_cov"), variable.name = "set", value.name = "n_sites")
dat.c[, set := factor(set, levels = c("Covered", "Evaluable"))]
dat.c[, `:=`(platform = plat_facet(dataset), cov = cov_facet(min_cov))]
drops.c <- thr_all[, .(lab = paste(sprintf("%s: -%.1f%% of sites", plat_short(dataset), pct_sites_dropped),
                                   collapse = "\n")), by = min_cov]
drops.c[, cov := cov_facet(min_cov)]
write_sd(dat.c, SRC_DIR, STEM, "c_sites")
write_sd(drops.c, SRC_DIR, STEM, "c_drops")
