#!/usr/bin/env Rscript
# gt_eval_sites.R  (dataset-agnostic GT evaluator, step 1/3)
#
# Collects the DISCOVERED A>G editing SITES (positions only, not levels) from a
# scRNAEditing Nextflow results dir, across every arm x method combination
# (caller x bam_type x strand_mode). Editing LEVELS are (re)quantified later by
# an independent pileup (gt_site_pileup.py) so the benchmark does not depend on
# any caller's own cov_ref/cov_alt or strand handling.
#
# Site source = pre100 tables by default (manuscript-final site set, before the
# 100%-editing removal step). By default only varType == "AG" rows are kept
# (A-to-I benchmark). Pass --all_vartypes true to keep all 12 mismatch types
# (for AG-enrichment among sites that pass the ground-truth criteria).
#
# Inputs (per arm): <results_dir>/<arm>/04_filtered/<caller>/<bam_type>/<strand_mode>/
#                   <arm>_<caller>_<bam_type>_<strand_mode>_<table>.tsv.gz
# Outputs (in --outdir):
#   gt_eval_discovery.tsv.gz   long: combo, caller, bam_type, strand_mode, arm,
#                                    chr, pos, strand, site [, varType]
#   gt_eval_union_sites.tsv    unique chr, pos, strand, site  (consensus strand)
#   gt_eval_union_sites.bed    chr, start0, end, site, ., strand  (for pileup)
#   gt_eval_site_vartype.tsv   (only with --all_vartypes) majority varType / site
#
# Usage:
#   Rscript gt_eval_sites.R --results_dir DIR --outdir DIR [--arms a,b,c] \
#           [--table pre100|filtered] [--all_vartypes true]

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || !length(x)) y else x
parse_cli <- function() {
  a <- commandArgs(trailingOnly = TRUE); out <- list()
  i <- 1L
  while (i <= length(a)) {
    if (grepl("^--", a[i])) { out[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
    else i <- i + 1L
  }
  out
}
args <- parse_cli()
results_dir <- args[["results_dir"]] %||% stop("--results_dir required")
outdir      <- args[["outdir"]]      %||% stop("--outdir required")
tbl         <- args[["table"]]       %||% "pre100"
all_vt      <- identical(tolower(args[["all_vartypes"]] %||% "false"), "true")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suffix <- if (tbl == "filtered") "filtered.txt.gz" else "pre100.tsv.gz"

# Arms: explicit or auto-detect (subdirs that contain 04_filtered/)
arms <- if (!is.null(args[["arms"]])) strsplit(args[["arms"]], ",")[[1]] else {
  d <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)
  d[vapply(d, function(x) dir.exists(file.path(results_dir, x, "04_filtered")), logical(1))]
}
if (!length(arms)) stop("No arms with 04_filtered/ under ", results_dir)
message("Arms: ", paste(arms, collapse = ", "))

read_gz <- function(f) {
  # R.utils may be absent; read via zcat.
  tryCatch(fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", header = TRUE),
           error = function(e) data.table())
}

recs <- list()
for (arm in arms) {
  files <- list.files(file.path(results_dir, arm, "04_filtered"),
                      pattern = paste0("_", suffix, "$"), recursive = TRUE,
                      full.names = TRUE)
  for (f in files) {
    # .../04_filtered/<caller>/<bam_type>/<strand_mode>[/snp_{off,on}]/<file>
    #
    # >>> snp_subtract IS A GRID DIMENSION AND MUST NOT BE POOLED. <<<
    # The publishDir layout gained a fourth `snp_{off,on}` level, but this
    # parser read only three, so snp_off and snp_on collapsed into a single
    # "combo" and every downstream count was the UNION of two different
    # methods. snp_subtract is one of the four declared GRID_DIMS in
    # publication/benchmark_grid/scripts/BENCHMARK_common.R, and the project
    # rule is that no unique combination may be silently dropped when
    # collating. Older runs without the level report snp_subtract = "na".
    parts <- strsplit(f, .Platform$file.sep)[[1]]
    k <- which(parts == "04_filtered")
    if (!length(k) || length(parts) < k + 3L) next
    caller <- parts[k + 1L]; bam_type <- parts[k + 2L]; strand_mode <- parts[k + 3L]
    snp_subtract <- if (length(parts) >= k + 4L && grepl("^snp_", parts[k + 4L]))
      sub("^snp_", "", parts[k + 4L]) else "na"
    d <- read_gz(f)
    if (!nrow(d)) next
    if (!all(c("chr", "pos", "varType") %in% names(d))) next
    if (!all_vt) d <- d[varType == "AG"]
    if (!nrow(d)) next
    if (!"strand" %in% names(d)) d[, strand := NA_character_]
    rec <- data.table(
      combo = paste(caller, bam_type, strand_mode, snp_subtract, sep = "__"),
      caller = caller, bam_type = bam_type, strand_mode = strand_mode,
      snp_subtract = snp_subtract,
      arm = arm, chr = as.character(d$chr), pos = as.integer(d$pos),
      strand = as.character(d$strand)
    )
    if (all_vt) rec[, varType := as.character(d$varType)]
    recs[[length(recs) + 1L]] <- rec
  }
}
if (!length(recs)) stop("No sites parsed from any ", suffix, " under ", results_dir,
                        " (are the downstream outputs populated?)")
disc <- rbindlist(recs, fill = TRUE)
disc[, site := paste0(chr, ":", pos)]
fwrite(disc, file.path(outdir, "gt_eval_discovery.tsv.gz"), sep = "\t")

# Consensus strand per site (majority; ties -> '+').
strand_vote <- disc[!is.na(strand) & strand %in% c("+", "-"),
                    .N, by = .(site, strand)]
consensus <- strand_vote[order(site, -N)][, .(strand = strand[1]), by = site]
union <- unique(disc[, .(site, chr, pos)])
union <- merge(union, consensus, by = "site", all.x = TRUE)
union[is.na(strand), strand := "+"]
setorder(union, chr, pos)
fwrite(union[, .(chr, pos, strand, site)],
       file.path(outdir, "gt_eval_union_sites.tsv"), sep = "\t")
# BED (0-based start)
fwrite(union[, .(chr, start = pos - 1L, end = pos, site, score = ".", strand)],
       file.path(outdir, "gt_eval_union_sites.bed"), sep = "\t", col.names = FALSE)

if (all_vt) {
  # Majority varType per site (ties broken by alphabetical order for stability)
  vt <- disc[!is.na(varType) & nzchar(varType), .N, by = .(site, varType)]
  setorder(vt, site, -N, varType)
  vt <- vt[, .SD[1], by = site]
  fwrite(vt[, .(site, varType, n_support = N)],
         file.path(outdir, "gt_eval_site_vartype.tsv"), sep = "\t")
  cat("Majority varType table:", nrow(vt), "sites\n")
  print(vt[, .N, by = varType][order(-N)])
}

cat("Discovered site rows:", nrow(disc),
    if (all_vt) " (all varTypes)" else " (AG only)", "\n")
cat("Method combinations:", uniqueN(disc$combo), "\n")
cat("Unique union sites:", nrow(union), "\n")
cat("Wrote:", file.path(outdir, "gt_eval_discovery.tsv.gz"), "\n")
cat("Wrote:", file.path(outdir, "gt_eval_union_sites.tsv"), "\n")
cat("Wrote:", file.path(outdir, "gt_eval_union_sites.bed"), "\n")
