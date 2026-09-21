#!/usr/bin/env Rscript
# requant_join.R — join strand-aware pileup re-quantification into a single
# per-site output file, ADDING canonical columns (originals preserved):
#     edit_level_pileup, cov_pileup, edited_reads, ref_reads
#
# Keyed on chr:pos. Sites absent from the pileup (zero strand-correct coverage)
# get NA. Idempotent: pre-existing requant columns are dropped and recomputed.
#
# Usage: requant_join.R --infile in.gz --pileup pileup.tsv --outfile out.gz
options(scipen = 999)
suppressPackageStartupMessages({ library(optparse); library(data.table) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--infile",  type = "character"),
  make_option("--pileup",  type = "character"),
  make_option("--outfile", type = "character")
)))

new_cols <- c("edit_level_pileup", "cov_pileup", "edited_reads", "ref_reads")

read_gz <- function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) return(data.table())
  tryCatch(fread(cmd = sprintf("zcat %s", shQuote(f)), sep = "\t", data.table = TRUE),
           error = function(e) data.table())
}

d <- read_gz(opt$infile)
if (!nrow(d) || !all(c("chr", "pos") %in% names(d))) {
  # nothing to annotate: copy input to output verbatim
  file.copy(opt$infile, opt$outfile, overwrite = TRUE)
  quit(status = 0)
}
d[, chr := as.character(chr)][, pos := as.integer(pos)]
for (nc in new_cols) if (nc %in% names(d)) d[, (nc) := NULL]
onames <- names(d)

pil <- if (file.exists(opt$pileup) && file.info(opt$pileup)$size > 0)
         fread(opt$pileup, sep = "\t", data.table = TRUE) else data.table()

if (nrow(pil)) {
  pil[, chr := as.character(chr)][, pos := as.integer(pos)]
  pil[, `:=`(edit_level_pileup = suppressWarnings(as.numeric(editing_level)),
             cov_pileup   = as.integer(coverage),
             edited_reads = as.integer(variant_cov),
             ref_reads    = as.integer(coverage) - as.integer(variant_cov))]
  look <- pil[, .(chr, pos, edit_level_pileup, cov_pileup, edited_reads, ref_reads)]
  setkey(look, chr, pos)
  d <- look[d, on = c("chr", "pos")]
} else {
  for (nc in new_cols) d[, (nc) := NA]
}
# Reorder to originals-then-canonical, but defensively: setcolorder requires a
# permutation of the CURRENT columns. Two hazards handled here:
#   (1) an already-requantified input leaves the new_cols in `onames`, so
#       c(onames, new_cols) lists them twice -> unique() collapses the dupes;
#   (2) a data.table join can rename collisions with an `i.` prefix -> keep only
#       names that actually exist, then append any remainder.
want <- unique(c(onames, new_cols))
want <- want[want %in% names(d)]
setcolorder(d, c(want, setdiff(names(d), want)))

# outfile ends in .gz -> fwrite gzips natively
fwrite(d, opt$outfile, sep = "\t", quote = FALSE, row.names = FALSE)
