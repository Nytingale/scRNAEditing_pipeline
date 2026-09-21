#!/usr/bin/env Rscript
# ============================================================================
# Per-cell downstream filtering
#
# Loads annotation GRanges ONCE, reads a collated per-cell variant TSV
# (with cell_id column), applies the same filtering logic as
# downstream_filter.R per cell, outputs a summary row per cell.
#
# Runs at both default and relaxed thresholds.
# ============================================================================
options(scipen = 999)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
})

option_list = list(
  make_option("--input",          type = "character", help = "Collated per-cell variants TSV.gz (with cell_id column)"),
  make_option("--caller",         type = "character", help = "Caller: bcftools, reditools3, gatk, scallele (scallele/reditools3/gatk share the single-alt-count 'ad' parser)"),
  make_option("--sample",         type = "character", help = "Sample identifier"),
  make_option("--genome_fasta",   type = "character", help = "Reference genome FASTA"),
  make_option("--snp_rdata",      type = "character", default = NULL, help = "Known SNPs GRanges RData"),
  make_option("--strand_annots",  type = "character", default = NULL, help = "Strand annotation RData"),
  make_option("--main_chroms",    type = "character", default = NULL, help = "Chromosomes to keep"),
  make_option("--min_cov",        type = "integer",   default = 10, help = "Default min coverage"),
  make_option("--min_alt",        type = "integer",   default = 3,  help = "Default min alt read count"),
  make_option("--min_edit_level", type = "double",    default = 0.01, help = "Default min editing level"),
  make_option("--output",         type = "character", help = "Output per-cell discovery summary TSV")
)
opt = parse_args(OptionParser(option_list = option_list))

BASES    = c("A", "C", "G", "T")
VARTYPES = as.character(outer(BASES, BASES, paste0))
VARTYPES = VARTYPES[!substr(VARTYPES, 1, 1) == substr(VARTYPES, 2, 2)]

log_msg = function(...) message("[percell_filter] ", ...)

# ── Load annotations once ────────────────────────────────────────────────────
log_msg("Loading annotations...")

snp_gr = NULL
if (!is.null(opt$snp_rdata) && file.exists(opt$snp_rdata)) {
  e_snp = new.env()
  load(opt$snp_rdata, envir = e_snp)
  if ("snp_gr" %in% ls(e_snp)) {
    snp_gr = e_snp$snp_gr
  } else {
    for (nm in ls(e_snp)) {
      obj = e_snp[[nm]]
      if (is.data.frame(obj) && "site" %in% colnames(obj)) {
        parts = data.table::tstrsplit(obj$site, ":", fixed = TRUE)
        snp_gr = GRanges(seqnames = parts[[1]], ranges = IRanges(start = as.integer(parts[[2]]), end = as.integer(parts[[2]])))
        snp_gr = sort(unique(snp_gr))
        break
      }
    }
  }
  rm(e_snp); gc()
  log_msg("  SNP GRanges: ", length(snp_gr), " positions")
}

plus_annot = NULL
minus_annot = NULL
if (!is.null(opt$strand_annots) && file.exists(opt$strand_annots)) {
  e_strand = new.env()
  load(opt$strand_annots, envir = e_strand)
  if ("plus.annot" %in% ls(e_strand))  plus_annot  = e_strand$plus.annot
  if ("minus.annot" %in% ls(e_strand)) minus_annot = e_strand$minus.annot
  rm(e_strand); gc()
  log_msg("  Strand annotations loaded: plus=", length(plus_annot), " minus=", length(minus_annot))
}

valid_chroms = NULL
if (!is.null(opt$main_chroms) && file.exists(opt$main_chroms)) {
  valid_chroms = readLines(opt$main_chroms)
  valid_chroms = valid_chroms[nchar(valid_chroms) > 0]
  log_msg("  Chromosome whitelist: ", length(valid_chroms), " chromosomes")
}

# ── Read collated per-cell variants ──────────────────────────────────────────
log_msg("Reading input: ", opt$input)
df = fread(opt$input, sep = "\t")
log_msg("  Total rows: ", nrow(df), " across ", length(unique(df$cell_id)), " cells")

if (nrow(df) == 0) {
  out = data.table(
    sample = opt$sample, cell_id = character(), caller = opt$caller,
    threshold_set = character(),
    total_sites = integer(), ag_sites = integer(), ag_pct = numeric(),
    total_after_filter = integer(), ag_after_filter = integer(), ag_pct_after_filter = numeric()
  )
  fwrite(out, opt$output, sep = "\t")
  log_msg("No variant data. Empty output written.")
  quit(status = 0)
}

# ── Standardize columns ─────────────────────────────────────────────────────
setnames(df, old = names(df), new = c("cell_id", "chr", "pos", "ref", "alt", "dp", "ad", "i16"))
df[, pos := as.integer(pos)]
df[, dp  := as.integer(dp)]

# Parse AD field (ref,alt format for bcftools)
if (opt$caller == "bcftools") {
  ad_parts = tstrsplit(df$ad, ",", fixed = TRUE, type.convert = TRUE)
  df[, cov_ref := ad_parts[[1]]]
  df[, cov_alt := ifelse(length(ad_parts) >= 2, ad_parts[[2]], 0L)]
} else {
  df[, cov_ref := dp - as.integer(gsub(",.*", "", ad))]
  df[, cov_alt := as.integer(gsub(",.*", "", ad))]
}

df[, varType := paste0(toupper(ref), toupper(alt))]
df[, site := paste0(chr, ":", pos)]

# ── Filter function (applied per cell) ───────────────────────────────────────
filter_cell = function(cell_df, min_cov, min_alt, min_edit_level) {
  # Main chroms
  if (!is.null(valid_chroms)) {
    cell_df = cell_df[chr %in% valid_chroms]
  }
  n_raw = nrow(cell_df)

  # Multi-mismatch removal
  cell_df = cell_df[varType %in% VARTYPES]

  # SNP removal
  if (!is.null(snp_gr) && nrow(cell_df) > 0) {
    gr = GRanges(seqnames = cell_df$chr, ranges = IRanges(start = cell_df$pos, end = cell_df$pos))
    hits = queryHits(findOverlaps(gr, snp_gr, type = "equal"))
    if (length(hits) > 0) cell_df = cell_df[-hits]
  }

  # Strand assignment (annotation-only for simplicity in per-cell mode)
  if (!is.null(plus_annot) && !is.null(minus_annot) && nrow(cell_df) > 0) {
    gr = GRanges(seqnames = cell_df$chr, ranges = IRanges(start = cell_df$pos, end = cell_df$pos))
    on_plus  = overlapsAny(gr, plus_annot)
    on_minus = overlapsAny(gr, minus_annot)
    cell_df[, strand := ifelse(on_plus & !on_minus, "+",
                        ifelse(!on_plus & on_minus, "-", NA_character_))]
    cell_df = cell_df[!is.na(strand)]

    # Flip to transcription direction
    comp = c(A = "T", T = "A", C = "G", G = "C")
    minus_rows = cell_df$strand == "-"
    if (any(minus_rows)) {
      cell_df[minus_rows, ref := comp[ref]]
      cell_df[minus_rows, alt := comp[alt]]
      cell_df[minus_rows, varType := paste0(ref, alt)]
    }
  }

  # Coverage filters
  if (nrow(cell_df) > 0) {
    cell_df[, cov := as.integer(cov_ref) + as.integer(cov_alt)]
    cell_df = cell_df[cov >= min_cov]
    cell_df = cell_df[!is.na(cov_alt) & cov_alt >= min_alt]
    cell_df[, edit_level := cov_alt / (cov_ref + cov_alt)]
    cell_df = cell_df[!is.na(edit_level) & edit_level >= min_edit_level]
    cell_df = cell_df[edit_level < 1]
  }

  # Compute metrics
  n_total = nrow(cell_df)
  n_ag    = sum(cell_df$varType == "AG", na.rm = TRUE)
  ag_pct  = if (n_total > 0) n_ag / n_total else NA_real_

  list(n_raw = n_raw, n_total = n_total, n_ag = n_ag, ag_pct = ag_pct)
}

# ── Process all cells ────────────────────────────────────────────────────────
cells = unique(df$cell_id)
log_msg("Processing ", length(cells), " cells at default + relaxed thresholds...")

results = list()
for (i in seq_along(cells)) {
  cid = cells[i]
  cell_df = df[cell_id == cid]

  if (i %% 500 == 0) log_msg("  Processed ", i, " / ", length(cells), " cells")

  # Default thresholds
  res_default = filter_cell(copy(cell_df), opt$min_cov, opt$min_alt, opt$min_edit_level)
  results[[length(results) + 1]] = data.table(
    sample = opt$sample, cell_id = cid, caller = opt$caller,
    threshold_set = "default",
    min_cov = opt$min_cov, min_alt = opt$min_alt, min_edit_level = opt$min_edit_level,
    raw_sites = res_default$n_raw,
    filtered_sites = res_default$n_total,
    ag_sites = res_default$n_ag,
    ag_pct = res_default$ag_pct
  )

  # Relaxed thresholds
  res_relaxed = filter_cell(copy(cell_df), 3L, 1L, 0.005)
  results[[length(results) + 1]] = data.table(
    sample = opt$sample, cell_id = cid, caller = opt$caller,
    threshold_set = "relaxed",
    min_cov = 3L, min_alt = 1L, min_edit_level = 0.005,
    raw_sites = res_relaxed$n_raw,
    filtered_sites = res_relaxed$n_total,
    ag_sites = res_relaxed$n_ag,
    ag_pct = res_relaxed$ag_pct
  )
}

out = rbindlist(results)
fwrite(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

n_default = out[threshold_set == "default"]
n_relaxed = out[threshold_set == "relaxed"]
log_msg("=== Summary (default thresholds) ===")
log_msg("  Cells with >0 filtered sites: ", sum(n_default$filtered_sites > 0), " / ", nrow(n_default))
log_msg("  Median filtered sites/cell:   ", median(n_default$filtered_sites))
log_msg("  Median AG sites/cell:         ", median(n_default$ag_sites))
log_msg("  Median AG purity:             ", round(median(n_default$ag_pct, na.rm = TRUE) * 100, 1), "%")
log_msg("=== Summary (relaxed thresholds) ===")
log_msg("  Cells with >0 filtered sites: ", sum(n_relaxed$filtered_sites > 0), " / ", nrow(n_relaxed))
log_msg("  Median filtered sites/cell:   ", median(n_relaxed$filtered_sites))
log_msg("  Median AG sites/cell:         ", median(n_relaxed$ag_sites))
log_msg("  Median AG purity:             ", round(median(n_relaxed$ag_pct, na.rm = TRUE) * 100, 1), "%")
log_msg("Output: ", opt$output)
log_msg("Done.")
