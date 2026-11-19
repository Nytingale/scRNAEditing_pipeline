#!/usr/bin/env Rscript
# requant_editing.R
# -----------------------------------------------------------------------------
# Post-hoc, method-independent re-quantification of editing level + coverage for
# every per-site output file under a results directory's 04_filtered/ tree.
#
# WHY: each caller reports edit_level = cov_alt/(cov_ref+cov_alt) from its own
# AD/base bookkeeping (different QC, BAQ, strand handling per method), so the
# published values are not strictly comparable and may be trusted verbatim by
# downstream users. We re-derive one canonical, strand-aware editing level and
# coverage directly from the *matching* bam_type BAM via the project's
# established stranded pileup parser (misc_scripts/parse_pileup_stranded.pl),
# and ADD them as new columns (originals preserved):
#     edit_level_pileup, cov_pileup, edited_reads, ref_reads
#
# The pileup editing level is a physical property of the reads at a site, so it
# is identical across strand_mode variants of the same (arm, bam_type) — only
# which sites *pass* differs by strand_mode. It depends on bam_type (dedup/recal
# change the reads), so we always use the BAM that produced each output.
#
# Efficiency: sites are pooled per (arm, bam_type); one mpileup pass per BAM.
#
# Usage:
#   Rscript requant_editing.R \
#     --results_dir human/pm/results_gt \
#     --fasta       human/ref_annots/hg38.analysisSet.fa \
#     [--libtype forward] [--arms NTC_ctrl,ADAR1_KO] [--bam_types sample,recal] \
#     [--pileup_sh human/scripts/gt_site_pileup.sh] [--min_bq 30 --min_mq 30] \
#     [--tmpdir /tmp/requant] [--backup] [--dry_run]
# -----------------------------------------------------------------------------
options(scipen = 999)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--results_dir", type = "character", help = "Dir containing <arm>/04_filtered/... trees"),
  make_option("--fasta",       type = "character", help = "Reference genome FASTA (indexed)"),
  make_option("--libtype",     type = "character", default = "forward", help = "forward|reverse|unstranded [forward]"),
  make_option("--arms",        type = "character", default = NULL, help = "Comma-list of arm subdirs (default: all with 04_filtered)"),
  make_option("--bam_types",   type = "character", default = NULL, help = "Comma-list of bam_labels to restrict to (default: all found)"),
  make_option("--pileup_sh",   type = "character",
              default = file.path(Sys.getenv("SCRNA_DIR", "."), "human/scripts/gt_site_pileup.sh"),
              help = "Path to gt_site_pileup.sh wrapper"),
  make_option("--min_bq",      type = "integer", default = 30L, help = "Min base quality [30]"),
  make_option("--min_mq",      type = "integer", default = 30L, help = "Min mapping quality [30]"),
  make_option("--max_depth",   type = "character", default = "1000000", help = "samtools mpileup -d [1e6]"),
  make_option("--tmpdir",      type = "character", default = NULL, help = "Scratch dir for BED + pileup (default: <results_dir>/.requant_tmp)"),
  make_option("--backup",      action = "store_true", default = FALSE, help = "Keep <file>.orig backups before rewriting"),
  make_option("--dry_run",     action = "store_true", default = FALSE, help = "List planned work; do not run pileup or rewrite")
)))

stopifnot(!is.null(opt$results_dir), !is.null(opt$fasta))
if (!dir.exists(opt$results_dir)) stop("results_dir not found: ", opt$results_dir)
if (!file.exists(opt$fasta))      stop("fasta not found: ", opt$fasta)
if (!file.exists(opt$pileup_sh))  stop("pileup_sh not found: ", opt$pileup_sh)

tmpdir <- if (is.null(opt$tmpdir)) file.path(opt$results_dir, ".requant_tmp") else opt$tmpdir
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── bam_type -> BAM path resolver (relative to an arm dir) ───────────────────
# Mirrors the pipeline's publishDir layout. Returns NA if no BAM is published.
resolve_bam <- function(arm_dir, arm, bam_label) {
  p <- function(...) file.path(arm_dir, ...)
  # dedup / recal variants share a physical BAM
  cand <- switch(bam_label,
    "sample"          = p("01_merged", paste0(arm, "_merged.bam")),
    "sampleRmDup"     = p("02_dedup",  paste0(arm, "_rmDup.bam")),
    "sampleUmiDup"    = p("02_dedup",  paste0(arm, "_umiDup.bam")),
    "sampleNoDup"     = p("02_dedup",  paste0(arm, "_noDup.bam")),
    "recal"           = p("02_gatk_preprocess", paste0(arm, "_recal.bam")),
    "recalq20"        = p("02_gatk_preprocess", paste0(arm, "_recal.bam")),
    "recaldefault"    = p("02_gatk_preprocess", paste0(arm, "_recal.bam")),
    "split"           = p("02_gatk_preprocess", paste0(arm, "_split.bam")),
    "goodCells"       = p("02_barcode_filtered", "goodCells",       paste0(arm, "_goodCells.bam")),
    "goodCellsRmDup"  = p("02_barcode_filtered", "goodCellsRmDup",  paste0(arm, "_goodCellsRmDup.bam")),
    "annotated"       = p("02_barcode_filtered", "annotated",       paste0(arm, "_annotated.bam")),
    "annotatedRmDup"  = p("02_barcode_filtered", "annotatedRmDup",  paste0(arm, "_annotatedRmDup.bam")),
    NA_character_)
  if (is.na(cand)) return(NA_character_)
  if (!file.exists(cand)) {
    # fallback: any BAM in the expected dir starting with arm_
    d <- dirname(cand)
    if (dir.exists(d)) {
      hits <- list.files(d, pattern = paste0("^", arm, ".*\\.bam$"), full.names = TRUE)
      if (length(hits) == 1) return(hits[1])
    }
    return(NA_character_)
  }
  cand
}

# ── enumerate per-site files under an arm's 04_filtered tree ─────────────────
# path parts after 04_filtered/: <caller>/<bam_label>/<strand_mode>[/snp_*]/<file>
list_site_files <- function(arm_dir) {
  root <- file.path(arm_dir, "04_filtered")
  if (!dir.exists(root)) return(data.table())
  files <- list.files(root, pattern = "_(filtered\\.txt|pre100\\.tsv)\\.gz$",
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(data.table())
  rel <- sub(paste0("^", root, "/"), "", files)
  parts <- strsplit(rel, "/", fixed = TRUE)
  caller    <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  bam_label <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  strand    <- vapply(parts, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))
  data.table(file = files, caller = caller, bam_label = bam_label, strand_mode = strand)
}

read_gz <- function(f) {
  if (file.info(f)$size == 0) return(data.table())
  tryCatch(fread(cmd = sprintf("zcat %s", shQuote(f)), sep = "\t", data.table = TRUE),
           error = function(e) { log_msg("  WARN read failed (", basename(f), "): ", e$message); data.table() })
}

# ── main ─────────────────────────────────────────────────────────────────────
arms <- if (!is.null(opt$arms)) strsplit(opt$arms, ",")[[1]] else basename(list.dirs(opt$results_dir, recursive = FALSE))
arms <- arms[dir.exists(file.path(opt$results_dir, arms, "04_filtered"))]
if (!length(arms)) stop("No arms with 04_filtered under ", opt$results_dir)
bt_filter <- if (!is.null(opt$bam_types)) strsplit(opt$bam_types, ",")[[1]] else NULL

log_msg("results_dir = ", opt$results_dir)
log_msg("arms        = ", paste(arms, collapse = ", "))
log_msg("libtype     = ", opt$libtype, "  MQ>=", opt$min_mq, " BQ>=", opt$min_bq,
        "  (mode: ", if (opt$dry_run) "DRY-RUN" else "APPLY", ")")

new_cols <- c("edit_level_pileup", "cov_pileup", "edited_reads", "ref_reads")
tot_files <- 0L; tot_updated <- 0L; skipped_nobam <- character()

for (arm in arms) {
  arm_dir <- file.path(opt$results_dir, arm)
  sf <- list_site_files(arm_dir)
  if (!nrow(sf)) { log_msg(arm, ": no per-site files"); next }
  if (!is.null(bt_filter)) sf <- sf[bam_label %in% bt_filter]
  if (!nrow(sf)) next
  log_msg("=== arm ", arm, " : ", nrow(sf), " files across bam_types {",
          paste(sort(unique(sf$bam_label)), collapse = ","), "} ===")

  for (bl in sort(unique(sf$bam_label))) {
    grp <- sf[bam_label == bl]
    bam <- resolve_bam(arm_dir, arm, bl)
    if (is.na(bam)) {
      log_msg("  [", bl, "] no published BAM -> skip ", nrow(grp), " files")
      skipped_nobam <- c(skipped_nobam, paste0(arm, "/", bl))
      next
    }

    # pooled union of sites across all files in this (arm, bam_type)
    site_tabs <- lapply(grp$file, read_gz)
    non_empty <- which(vapply(site_tabs, nrow, integer(1)) > 0)
    if (!length(non_empty)) { log_msg("  [", bl, "] all files empty -> skip"); next }
    keys <- unique(rbindlist(lapply(site_tabs[non_empty], function(d) {
      if (!all(c("chr", "pos") %in% names(d))) return(data.table())
      unique(d[, .(chr = as.character(chr), pos = as.integer(pos))])
    }), fill = TRUE))
    if (!nrow(keys) || !all(c("chr", "pos") %in% names(keys))) { log_msg("  [", bl, "] no chr/pos keys -> skip"); next }
    keys <- keys[!is.na(chr) & !is.na(pos)]
    if (!nrow(keys)) { log_msg("  [", bl, "] no chr/pos keys -> skip"); next }
    log_msg("  [", bl, "] BAM=", basename(bam), "  union sites=", nrow(keys),
            "  files=", nrow(grp))

    if (opt$dry_run) { tot_files <- tot_files + nrow(grp); next }

    # BED (0-based half-open) + single mpileup pass
    bed <- file.path(tmpdir, sprintf("%s_%s.bed", arm, bl))
    out <- file.path(tmpdir, sprintf("%s_%s.pileup.tsv", arm, bl))
    b <- keys[order(chr, pos)][, .(chr, start = pos - 1L, end = pos)]
    fwrite(b, bed, sep = "\t", col.names = FALSE)
    cmd <- sprintf("bash %s --bam %s --bed %s --fasta %s --out %s --libtype %s --min-bq %d --min-mq %d --max-depth %s",
                   shQuote(opt$pileup_sh), shQuote(bam), shQuote(bed), shQuote(opt$fasta),
                   shQuote(out), opt$libtype, opt$min_bq, opt$min_mq, opt$max_depth)
    rc <- system(cmd)
    if (rc != 0 || !file.exists(out)) { log_msg("  [", bl, "] pileup FAILED (rc=", rc, ") -> skip"); next }

    pil <- fread(out, sep = "\t", data.table = TRUE)
    if (!nrow(pil)) { log_msg("  [", bl, "] pileup empty -> skip"); next }
    # perl cols: chr pos ref coverage variant_cov A C G T a c g t ref_count_plus ref_count_minus editing_level
    pil[, chr := as.character(chr)][, pos := as.integer(pos)]
    pil[, `:=`(edit_level_pileup = suppressWarnings(as.numeric(editing_level)),
               cov_pileup   = as.integer(coverage),
               edited_reads = as.integer(variant_cov),
               ref_reads    = as.integer(coverage) - as.integer(variant_cov))]
    look <- pil[, .(chr, pos, edit_level_pileup, cov_pileup, edited_reads, ref_reads)]
    setkey(look, chr, pos)

    # join back into each file, add columns, rewrite
    for (i in seq_len(nrow(grp))) {
      f <- grp$file[i]
      d <- site_tabs[[ which(grp$file == f)[1] ]]
      if (!nrow(d) || !all(c("chr", "pos") %in% names(d))) { tot_files <- tot_files + 1L; next }
      d[, chr := as.character(chr)][, pos := as.integer(pos)]
      for (nc in new_cols) if (nc %in% names(d)) d[, (nc) := NULL]  # idempotent re-run
      onames <- names(d)                  # original column order (pre-join)
      d <- look[d, on = c("chr", "pos")]  # right join: keep all site rows
      setcolorder(d, c(onames, new_cols)) # original order, requant cols appended
      if (opt$backup && !file.exists(paste0(f, ".orig"))) file.copy(f, paste0(f, ".orig"))
      tmpf <- paste0(f, ".rq.tmp.gz")  # must end in .gz so fwrite gzips
      fwrite(d, tmpf, sep = "\t", quote = FALSE, row.names = FALSE)
      file.rename(tmpf, f)
      tot_files <- tot_files + 1L; tot_updated <- tot_updated + 1L
    }
    log_msg("  [", bl, "] updated ", nrow(grp), " files")
  }
}

log_msg("DONE. files seen=", tot_files, "  updated=", tot_updated,
        if (length(skipped_nobam)) paste0("  no-BAM skipped: ", paste(skipped_nobam, collapse = ", ")) else "")
