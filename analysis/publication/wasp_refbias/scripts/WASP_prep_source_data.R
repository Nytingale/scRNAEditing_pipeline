#!/usr/bin/env Rscript
# ============================================================================
# WASP_prep_source_data.R — source-data generation for the WASP reference-bias
# supplement (Reviewer 2, Major 4). NO PLOTTING HAPPENS HERE.
# ============================================================================
# Stage 2 of 2 in the source-data chain:
#   Stage 1  WASP_prep_pileups.sh   FASTQ/BAM -> merged BAMs -> stranded
#                                   pileups at REDIportal A-to-I sites
#                                   (+ global read-loss counts)
#   Stage 2  THIS SCRIPT            pileups -> paired per-site table (the
#                                   shared plot-ready input) + report tables
#
# The comparison: strand-aware editing levels at REDIportal A-to-I sites in
#   ALL reads          (standard STAR alignment)        = biased estimate
#   WASP-passed reads  (vW==1 or no vW tag)             = bias-corrected
# STAR's WASP flips the allele at each variant, re-maps, and flags reads whose
# mapping is allele-dependent (vW>1) — precisely the reference-biased reads.
# If reference bias suppresses editing, removing them should RAISE editing
# levels, because reference-supporting reads are preferentially removed.
#
# Usage:
#   Rscript WASP_prep_source_data.R [DATASET] [MINCOV]
# Defaults: DATASET = jzlxv3, MINCOV = the registry value for that dataset.
# Every path is resolved from scripts/WASP_datasets.tsv, so adding a dataset
# is a row in that file rather than an edit here.
#
# Reads (from the dataset's analysis dir, all produced by Stage 1):
#   pileup_all.tsv, pileup_pass.tsv, readloss.tsv
#   01_starsolobams_wasp_filtered/wasp_filter_summary.csv  (per-unit vW)
#
# Writes, all under a <dataset> subfolder:
#   data/<ds>/WASP_persite_evaluable.tsv.gz  paired evaluable sites (figure input)
#   report/<ds>/WASP_by_threshold.tsv        coverage-threshold sensitivity sweep
#   report/<ds>/WASP_removed_alleles.tsv     allele composition of removed reads
#   report/<ds>/WASP_readloss.tsv            read loss, global AND WASP-evaluated
#   report/<ds>/WASP_summary.tsv             headline metrics (long format)
#   source_data/<ds>/SOURCE_DATA_MANIFEST.tsv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

this_dir <- local({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
})
DELIV_DIR <- normalizePath(file.path(this_dir, ".."), mustWork = FALSE)
REPO_DIR  <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else "jzlxv3"

reg <- fread(file.path(this_dir, "WASP_datasets.tsv"))
cfg <- reg[dataset == DATASET]
if (nrow(cfg) != 1L) {
  stop("Unknown dataset '", DATASET, "'. Known: ",
       paste(reg$dataset, collapse = ", "))
}
MINCOV <- if (length(args) >= 2) as.integer(args[2]) else as.integer(cfg$mincov)

ANALYSIS_DIR <- file.path(REPO_DIR, "human", DATASET, "wasp_refbias")
DATA_DIR   <- file.path(DELIV_DIR, "data",        DATASET)
REPORT_DIR <- file.path(DELIV_DIR, "report",      DATASET)
SRC_DIR    <- file.path(DELIV_DIR, "source_data", DATASET)
for (d in c(DATA_DIR, REPORT_DIR, SRC_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

f_all  <- file.path(ANALYSIS_DIR, "pileup_all.tsv")
f_pass <- file.path(ANALYSIS_DIR, "pileup_pass.tsv")
f_rl   <- file.path(ANALYSIS_DIR, "readloss.tsv")
f_vw   <- file.path(REPO_DIR, "human", DATASET,
                    "01_starsolobams_wasp_filtered", "wasp_filter_summary.csv")
stopifnot(file.exists(f_all), file.exists(f_pass))

cat("=== WASP reference-bias source data ===\n")
cat("dataset      :", DATASET, "(", cfg$label, ")\n")
cat("analysis dir :", ANALYSIS_DIR, "\n")
cat("min coverage :", MINCOV, "\n\n")

# ── Per-site pileups ────────────────────────────────────────────────────────
# Strand-aware columns from the project pileup engine (gt_site_pileup.sh):
#   coverage    = ref + edited reads on the transcribed strand at the site
#   variant_cov = edited-allele (G on +, C on -) read count
# Only these columns are read; the 12 per-base columns are not needed and
# skipping them keeps the merge well inside memory.
read_pileup <- function(f) {
  d <- fread(f, sep = "\t", header = TRUE,
             select = c("chr", "pos", "ref", "coverage", "variant_cov",
                        "editing_level"))
  d[, editing_level := suppressWarnings(as.numeric(editing_level))]
  d[, ref_at_site := coverage - variant_cov]
  d[]
}

cat("--- reading pileups ---\n")
da <- read_pileup(f_all)
dp <- read_pileup(f_pass)
cat("  all  reads :", format(nrow(da), big.mark = ","), "sites\n")
cat("  WASP-passed:", format(nrow(dp), big.mark = ","), "sites\n")

setnames(da, c("coverage", "variant_cov", "ref_at_site", "editing_level"),
             c("cov_all", "alt_all", "ref_all", "EL_all"))
setnames(dp, c("coverage", "variant_cov", "ref_at_site", "editing_level"),
             c("cov_pass", "alt_pass", "ref_pass", "EL_pass"))

# Base set = sites adequately covered WITHOUT WASP. Left-join the WASP-passed
# pileup; a site absent there had ALL its reads removed, hence cov_pass = 0.
m <- merge(da, dp[, .(chr, pos, cov_pass, alt_pass, ref_pass, EL_pass)],
           by = c("chr", "pos"), all.x = TRUE)
for (col in c("cov_pass", "alt_pass", "ref_pass")) {
  m[is.na(get(col)), (col) := 0L]
}

# ── Coverage thresholds ─────────────────────────────────────────────────────
# MINCOV is the PRIMARY threshold: the one the headline numbers quote. The
# figures show all three side by side, because the answer's dependence on the
# floor is itself a result — the removed-read gap collapses from 1.08 pp at
# cov >= 10 to 0.03 pp (p = 0.24) at cov >= 50.
thresholds <- sort(unique(c(10L, 20L, 50L, MINCOV)))
EV_FLOOR   <- min(thresholds)

# ── Shared plot-ready input: the paired evaluable set ───────────────────────
# "Evaluable" = coverage >= floor in BOTH arms, so editing can be compared
# fairly, and a finite editing level on both sides.
#
# WRITTEN AT THE LOWEST THRESHOLD, not at MINCOV. The thresholds nest — the
# cov >= 10 evaluable set contains the cov >= 20 set, which contains cov >= 50 —
# so one table plus the cov_all / cov_pass columns lets the figures facet by
# threshold with a filter, instead of storing three overlapping copies that
# could drift apart.
ev <- m[cov_all >= EV_FLOOR & cov_pass >= EV_FLOOR &
          is.finite(EL_all) & is.finite(EL_pass),
        .(chr, pos, ref, cov_all, alt_all, ref_all, EL_all,
          cov_pass, alt_pass, ref_pass, EL_pass)]
ev[, dEL := EL_pass - EL_all]
setorder(ev, chr, pos)
fwrite(ev, file.path(DATA_DIR, "WASP_persite_evaluable.tsv.gz"), sep = "\t")
cat("  evaluable  :", format(nrow(ev), big.mark = ","),
    "paired sites at cov >=", EV_FLOOR, "(filter up for higher floors)\n\n")

# ── Coverage-threshold sensitivity sweep ────────────────────────────────────
# One pass per threshold produces BOTH the per-site statistics and the
# removed-read allele composition. The allele rates used to be computed at
# MINCOV alone, so the threshold dependence of the bias signature — the thing
# the figures now show — could not be read off any table.
cat("--- threshold sweep ---\n")
sweep <- lapply(thresholds, function(th) {
  base <- m[cov_all >= th]
  e <- base[cov_pass >= th & is.finite(EL_all) & is.finite(EL_pass)]
  d <- e$EL_pass - e$EL_all

  rem_ref <- sum(pmax(base$ref_all - base$ref_pass, 0))
  rem_alt <- sum(pmax(base$alt_all - base$alt_pass, 0))
  tot_ref <- sum(base$ref_all)
  tot_alt <- sum(base$alt_all)
  rate_ref <- if (tot_ref > 0) rem_ref / tot_ref else NA_real_
  rate_alt <- if (tot_alt > 0) rem_alt / tot_alt else NA_real_

  # Two-proportion test on the gap. With hundreds of millions of reference
  # reads any gap is "significant" at low thresholds, so the CI is what to
  # read, not the p-value; at cov >= 50 the CI straddles zero.
  ht <- tryCatch(prop.test(c(rem_ref, rem_alt), c(tot_ref, tot_alt)),
                 error = function(x) NULL)

  list(
    thr = data.table(
      min_cov           = th,
      n_sites_covered   = nrow(base),
      n_sites_evaluable = nrow(e),
      n_dropped_by_wasp = sum(base$cov_pass < th),
      pct_sites_dropped = round(100 * mean(base$cov_pass < th), 2),
      median_EL_all     = round(median(e$EL_all), 4),
      median_EL_wasp    = round(median(e$EL_pass), 4),
      mean_EL_all       = round(mean(e$EL_all), 4),
      mean_EL_wasp      = round(mean(e$EL_pass), 4),
      median_dEL        = round(median(d), 4),
      mean_dEL          = round(mean(d), 4),
      pct_sites_EL_up   = round(100 * mean(d > 0), 2),
      pct_sites_EL_down = round(100 * mean(d < 0), 2),
      removal_rate_ref  = round(rate_ref, 6),
      removal_rate_alt  = round(rate_alt, 6),
      gap_pp            = round(100 * (rate_ref - rate_alt), 3),
      gap_ci_lo_pp      = if (is.null(ht)) NA_real_ else round(100 * ht$conf.int[1], 3),
      gap_ci_hi_pp      = if (is.null(ht)) NA_real_ else round(100 * ht$conf.int[2], 3),
      gap_p             = if (is.null(ht)) NA_real_ else ht$p.value,
      wilcoxon_p        = tryCatch(
        wilcox.test(e$EL_pass, e$EL_all, paired = TRUE)$p.value,
        error = function(x) NA_real_)
    ),
    alleles = data.table(
      min_cov       = th,
      allele        = c("reference (A)", "edited (G)"),
      reads_all     = c(tot_ref, tot_alt),
      reads_removed = c(rem_ref, rem_alt),
      reads_kept    = c(tot_ref - rem_ref, tot_alt - rem_alt),
      removal_rate  = c(rate_ref, rate_alt)
    )
  )
})

by_thr <- rbindlist(lapply(sweep, `[[`, "thr"))
fwrite(by_thr, file.path(REPORT_DIR, "WASP_by_threshold.tsv"), sep = "\t")
print(by_thr[, .(min_cov, n_sites_covered, n_sites_evaluable, pct_sites_dropped,
                 mean_dEL, gap_pp, gap_ci_lo_pp, gap_ci_hi_pp)])

# ── Allele composition of WASP-removed reads (the direct bias signature) ────
# Under reference bias, WASP should remove reference-allele reads at a HIGHER
# rate than edited-allele reads at editing sites. Long format, one block per
# threshold, so the figure can facet on min_cov.
removed <- rbindlist(lapply(sweep, `[[`, "alleles"))
fwrite(removed, file.path(REPORT_DIR, "WASP_removed_alleles.tsv"), sep = "\t")
cat("\n--- removed-read allele composition, by threshold ---\n")
print(removed)

base <- m[cov_all >= MINCOV]   # primary threshold, for the read-loss table below

# ── Read loss: global AND among WASP-EVALUATED reads ───────────────────────
# IMPORTANT. The headline "~14.8% of reads removed" is a GENOME-WIDE average
# and is heavily diluted: STAR's WASP only evaluates reads that overlap a VCF
# variant, and reads overlapping nothing pass untouched. The cost that actually
# applies at editing sites is the loss among EVALUATED reads, which is far
# larger. Both are reported so the supplement cannot be read as understating
# the coverage cost.
rl <- if (file.exists(f_rl)) fread(f_rl) else NULL
rows <- list()
if (!is.null(rl)) {
  rows <- c(rows, list(data.table(
    scope  = "global (all aligned reads)",
    metric = rl$metric,
    value  = as.character(rl$value))))
}
if (file.exists(f_vw)) {
  vw <- fread(f_vw)
  tot_reads <- sum(vw$total_reads)
  novar     <- sum(vw$no_variant_overlap)
  passed    <- sum(vw$wasp_passed)
  failed    <- sum(vw$wasp_failed)
  evaluated <- passed + failed
  rows <- c(rows, list(data.table(
    scope  = "WASP-evaluated reads (overlap >=1 REDIportal variant)",
    metric = c("total_reads", "no_variant_overlap", "reads_evaluated",
               "reads_passed", "reads_failed", "frac_failed_of_evaluated",
               "frac_failed_of_total"),
    value  = c(as.character(tot_reads), as.character(novar),
               as.character(evaluated), as.character(passed),
               as.character(failed),
               sprintf("%.6f", failed / evaluated),
               sprintf("%.6f", failed / tot_reads)))))
}
# Coverage actually lost at the REDIportal sites being analysed.
rows <- c(rows, list(data.table(
  scope  = paste0("coverage at REDIportal A-to-I sites (cov >= ", MINCOV, ")"),
  metric = c("coverage_all", "coverage_wasp_pass", "frac_coverage_lost"),
  value  = c(as.character(sum(base$cov_all)),
             as.character(sum(base$cov_pass)),
             sprintf("%.6f", 1 - sum(base$cov_pass) / sum(base$cov_all))))))
readloss <- rbindlist(rows, fill = TRUE)
fwrite(readloss, file.path(REPORT_DIR, "WASP_readloss.tsv"), sep = "\t")
cat("\n--- read / coverage loss ---\n")
print(readloss)

# ── Headline summary ────────────────────────────────────────────────────────
prim <- by_thr[min_cov == MINCOV]
summary_long <- rbindlist(list(
  data.table(metric = "dataset",                  value = DATASET),
  data.table(metric = "dataset_label",            value = as.character(cfg$label)),
  data.table(metric = "platform",                 value = as.character(cfg$platform)),
  data.table(metric = "pileup_libtype",           value = as.character(cfg$libtype)),
  # min_coverage is the contract with the figures: WASP_common.R reads this row
  # back rather than hardcoding a threshold, so a panel can never state a
  # coverage floor different from the one its own data was built with.
  data.table(metric = "min_coverage",             value = as.character(MINCOV)),
  data.table(metric = "n_sites_evaluable",        value = as.character(prim$n_sites_evaluable)),
  data.table(metric = "median_EL_no_wasp",        value = as.character(prim$median_EL_all)),
  data.table(metric = "median_EL_wasp",           value = as.character(prim$median_EL_wasp)),
  data.table(metric = "mean_EL_no_wasp",          value = as.character(prim$mean_EL_all)),
  data.table(metric = "mean_EL_wasp",             value = as.character(prim$mean_EL_wasp)),
  data.table(metric = "median_delta_EL",          value = as.character(prim$median_dEL)),
  data.table(metric = "mean_delta_EL",            value = as.character(prim$mean_dEL)),
  data.table(metric = "pct_sites_EL_increased",   value = as.character(prim$pct_sites_EL_up)),
  data.table(metric = "pct_sites_EL_decreased",   value = as.character(prim$pct_sites_EL_down)),
  data.table(metric = "wilcoxon_paired_p",        value = format(prim$wilcoxon_p, scientific = TRUE)),
  # Index by threshold, not by row position: `removed` is now long over all
  # thresholds, so removed$removal_rate[1] would silently quote cov >= 10.
  data.table(metric = "ref_reads_removed_frac",
             value = sprintf("%.4f", removed[min_cov == MINCOV &
                                               allele == "reference (A)", removal_rate])),
  data.table(metric = "edited_reads_removed_frac",
             value = sprintf("%.4f", removed[min_cov == MINCOV &
                                               allele == "edited (G)", removal_rate])),
  data.table(metric = "removal_gap_pp",           value = sprintf("%.3f", prim$gap_pp)),
  data.table(metric = "removal_gap_ci_pp",
             value = sprintf("%.3f to %.3f", prim$gap_ci_lo_pp, prim$gap_ci_hi_pp)),
  data.table(metric = "evaluable_set_floor",      value = as.character(EV_FLOOR)),
  data.table(metric = "thresholds_swept",         value = paste(thresholds, collapse = ",")),
  readloss[, .(metric = paste0(scope, " | ", metric), value)]
), fill = TRUE)
fwrite(summary_long, file.path(REPORT_DIR, "WASP_summary.tsv"), sep = "\t")

# ── Source-data manifest ────────────────────────────────────────────────────
manifest <- data.table(
  file = c("pileup_all.tsv", "pileup_pass.tsv", "readloss.tsv",
           "wasp_filter_summary.csv", "WASP_persite_evaluable.tsv.gz"),
  role = c("stranded pileup at REDIportal A-to-I sites, ALL reads (biased estimate)",
           "stranded pileup at REDIportal A-to-I sites, WASP-passed reads only",
           "global aligned-read counts before/after WASP filtering",
           "per-unit vW tag breakdown (evaluated / passed / failed)",
           "paired evaluable per-site table; shared input to all WASP figures"),
  stage = c("WASP_prep_pileups.sh", "WASP_prep_pileups.sh",
            "WASP_prep_pileups.sh", "WASP_filter_bams.sh",
            "WASP_prep_source_data.R"),
  origin_dir = c(ANALYSIS_DIR, ANALYSIS_DIR, ANALYSIS_DIR,
                 dirname(f_vw), DATA_DIR)
)
manifest[, dataset := DATASET]
manifest[, exists := file.exists(file.path(origin_dir, file))]
setcolorder(manifest, c("dataset", "file", "role", "stage", "origin_dir", "exists"))
fwrite(manifest, file.path(SRC_DIR, "SOURCE_DATA_MANIFEST.tsv"), sep = "\t")

cat("\nWrote source data + report tables for", nrow(ev), "evaluable sites.\n")
cat("Next: bash scripts/run_WASP_figures.sh", DATASET, "\n")
