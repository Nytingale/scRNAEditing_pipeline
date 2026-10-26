#!/usr/bin/env Rscript
# =============================================================================
# audit_dedup_editing_relevance.R
#
# PI concern: is UMI/coordinate dedup needed, or are duplicates negligible for
# editing calculations?
#
# Two layers:
#   1. GLOBAL — MarkDuplicates PERCENT_DUPLICATION + UMI retention
#      (sample BAM read count vs umiDup / rmDup).
#   2. SITE-LEVEL — at sites that enter editing discovery (union of AG sites
#      called on the Unprocessed `sample` bam_type, or the v3 truth TRUE set),
#      compare samtools depth on merged vs rmDup BAMs:
#         frac_dup_depth = 1 - depth_rmDup / depth_sample
#      Restricted to positions with depth_sample >= min_cov.
#
# Writes under --outdir (default publication/teratoma_gt/report/pm_gt/):
#   DEDUP_global_rates.tsv
#   DEDUP_site_depth_frac.tsv          per-arm summary
#   DEDUP_site_depth_frac_detail.tsv.gz  optional per-site (if --write_sites)
#   DEDUP_f1_delta_by_bamtype.tsv      from gt_tier_scoring
# =============================================================================

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
parse_cli <- function() {
  a <- commandArgs(trailingOnly = TRUE); out <- list(); i <- 1L
  while (i <= length(a)) {
    if (grepl("^--", a[i])) { out[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L
    } else i <- i + 1L
  }
  out
}
args <- parse_cli()

BASE <- args[["base"]] %||%
  Sys.getenv("SCRNA_DIR", ".")
outdir <- args[["outdir"]] %||%
  file.path(BASE, "publication/teratoma_gt/report/pm_gt")
min_cov <- as.integer(args[["min_cov"]] %||% 10L)
write_sites <- identical(args[["write_sites"]] %||% "false", "true")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# ── 1. Global rates ────────────────────────────────────────────────────────
parse_picard <- function(f, dataset, arm) {
  lines <- readLines(f)
  i <- which(grepl("^LIBRARY", lines))[1]
  if (is.na(i)) return(NULL)
  hdr <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
  val <- strsplit(lines[i + 1], "\t", fixed = TRUE)[[1]]
  # Picard occasionally emits a trailing empty field or omits ESTIMATED_LIBRARY_SIZE;
  # align lengths before naming.
  if (length(val) < length(hdr)) val <- c(val, rep(NA_character_, length(hdr) - length(val)))
  if (length(val) > length(hdr)) val <- val[seq_along(hdr)]
  names(val) <- hdr
  data.table(
    dataset = dataset, arm = arm,
    unpaired_reads = as.numeric(val[["UNPAIRED_READS_EXAMINED"]]),
    read_pairs = as.numeric(val[["READ_PAIRS_EXAMINED"]]),
    unpaired_dups = as.numeric(val[["UNPAIRED_READ_DUPLICATES"]]),
    pair_dups = as.numeric(val[["READ_PAIR_DUPLICATES"]]),
    percent_duplication = as.numeric(val[["PERCENT_DUPLICATION"]]))
}

global <- list()
for (spec in list(
  list(ds = "teratoma_merged", root = file.path(BASE, "human/pm/results_gt_v2")),
  list(ds = "replogle",        root = file.path(BASE, "human/replogle/results_gt"))
)) {
  mets <- Sys.glob(file.path(spec$root, "*/02_dedup/*_dup_metrics.txt"))
  for (f in mets) {
    arm <- basename(dirname(dirname(f)))
    global[[length(global) + 1L]] <- parse_picard(f, spec$ds, arm)
    # UMI retention via idxstats if BAMs exist
    umi <- file.path(dirname(f), paste0(arm, "_umiDup.bam"))
    rmd <- file.path(dirname(f), paste0(arm, "_rmDup.bam"))
    merged <- file.path(spec$root, arm, "01_merged", paste0(arm, "_merged.bam"))
    if (!file.exists(merged))
      merged <- file.path(spec$root, arm, "00_validated",
                          paste0("merged_", arm, "_uniq.bam"))
    cnt <- function(b) {
      if (!file.exists(b)) return(NA_real_)
      as.numeric(system(sprintf("samtools idxstats %s | awk '{s+=$3} END{print s}'",
                                shQuote(b)), intern = TRUE)[1])
    }
    g <- global[[length(global)]]
    g[, `:=`(n_merged = cnt(merged), n_rmDup = cnt(rmd), n_umiDup = cnt(umi))]
    g[, `:=`(
      frac_removed_coord = 1 - n_rmDup / n_merged,
      frac_removed_umi   = 1 - n_umiDup / n_merged)]
    global[[length(global)]] <- g
  }
}
glob_dt <- rbindlist(global, fill = TRUE)
fwrite(glob_dt, file.path(outdir, "DEDUP_global_rates.tsv"), sep = "\t")
msg("global rates: ", nrow(glob_dt), " arms")

# ── 2. F1 deltas from scoring ──────────────────────────────────────────────
score_f <- file.path(BASE, "human/pm/gt_eval/true_def_v3_full/gt_tier_scoring.tsv")
if (file.exists(score_f)) {
  sc <- fread(score_f)
  sc[, f1 := fifelse(precision + recall_vs_union > 0,
                     2 * precision * recall_vs_union / (precision + recall_vs_union),
                     NA_real_)]
  wide <- dcast(sc[bam_type %in% c("sample", "sampleRmDup", "sampleUmiDup")],
                tier + caller + strand_mode + snp_subtract ~ bam_type,
                value.var = "f1")
  if (all(c("sample", "sampleRmDup") %in% names(wide))) {
    wide[, `:=`(
      delta_sample_minus_RmDup = sample - sampleRmDup,
      delta_sample_minus_Umi   = sample - sampleUmiDup)]
  }
  fwrite(wide, file.path(outdir, "DEDUP_f1_delta_by_bamtype.tsv"), sep = "\t")
  msg("F1 deltas written")
}

# ── 3. Site-level depth fraction ───────────────────────────────────────────
# Site list: TRUE sites from v3 truth (editing-relevant by definition).
truth_f <- file.path(BASE, "human/pm/gt_eval/true_def_v3_full/gt_truth_sites.tsv.gz")
if (!file.exists(truth_f)) {
  msg("no truth set; skipping site-level depth")
  quit(status = 0)
}
truth <- fread(cmd = paste("gzip -cd", shQuote(truth_f)), sep = "\t",
               showProgress = FALSE)
sites <- truth[class == "TRUE", .(chr = sub(":.*", "", site),
                                  pos = as.integer(sub(".*:", "", site)),
                                  site)]
bed <- tempfile(fileext = ".bed")
fwrite(unique(sites)[, .(chr, start = pos - 1L, end = pos)], bed,
       sep = "\t", col.names = FALSE)

depth_one <- function(bam, bed) {
  if (!file.exists(bam)) return(NULL)
  cmd <- sprintf("samtools depth -a -b %s %s", shQuote(bed), shQuote(bam))
  raw <- fread(cmd = cmd, col.names = c("chr", "pos", "depth"),
               showProgress = FALSE)
  raw
}

summaries <- list()
details <- list()
root <- file.path(BASE, "human/pm/results_gt_v2")
for (arm in c("AAVS1_ctrl", "ADAR1_KO", "ADARB1_KO", "ADARB2_KO", "NTC_ctrl")) {
  merged <- file.path(root, arm, "01_merged", paste0(arm, "_merged.bam"))
  rmd    <- file.path(root, arm, "02_dedup",  paste0(arm, "_rmDup.bam"))
  if (!file.exists(merged) || !file.exists(rmd)) {
    msg("skip ", arm, ": missing BAM")
    next
  }
  msg("depth ", arm)
  d0 <- depth_one(merged, bed)
  d1 <- depth_one(rmd, bed)
  if (is.null(d0) || is.null(d1)) next
  m <- merge(d0, d1, by = c("chr", "pos"), suffixes = c("_sample", "_rmDup"))
  m <- m[depth_sample >= min_cov]
  m[, frac_dup_depth := fifelse(depth_sample > 0,
                                1 - depth_rmDup / depth_sample, NA_real_)]
  summaries[[length(summaries) + 1L]] <- m[, .(
    dataset = "teratoma_merged", arm = arm, n_sites = .N,
    median_frac_dup_depth = median(frac_dup_depth, na.rm = TRUE),
    mean_frac_dup_depth   = mean(frac_dup_depth, na.rm = TRUE),
    q25 = quantile(frac_dup_depth, 0.25, na.rm = TRUE),
    q75 = quantile(frac_dup_depth, 0.75, na.rm = TRUE),
    median_depth_sample = median(depth_sample),
    median_depth_rmDup  = median(depth_rmDup))]
  if (write_sites) {
    m[, arm := arm]
    details[[length(details) + 1L]] <- m
  }
}
if (length(summaries)) {
  sum_dt <- rbindlist(summaries)
  fwrite(sum_dt, file.path(outdir, "DEDUP_site_depth_frac.tsv"), sep = "\t")
  msg("site-level summary: ", nrow(sum_dt), " arms; median frac_dup_depth range ",
      round(min(sum_dt$median_frac_dup_depth), 3), "-",
      round(max(sum_dt$median_frac_dup_depth), 3))
}
if (write_sites && length(details)) {
  det <- rbindlist(details, fill = TRUE)
  fwrite(det, file.path(outdir, "DEDUP_site_depth_frac_detail.tsv.gz"), sep = "\t")
}
msg("done -> ", outdir)
