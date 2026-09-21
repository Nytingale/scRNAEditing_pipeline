#!/usr/bin/env Rscript
# requant_sites_to_bed.R — union of (chr,pos) across one or more per-site files
# into a 0-based BED for samtools mpileup -l. Robust to differing column sets
# (uses named 'chr'/'pos' columns). Writes an empty file if no sites.
#
# Usage: requant_sites_to_bed.R --files a.gz,b.gz --bed sites.bed
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--files", type = "character"),
  make_option("--bed",   type = "character")
)))
read_gz <- function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) return(data.table())
  tryCatch(fread(cmd = sprintf("zcat %s", shQuote(f)), sep = "\t", data.table = TRUE),
           error = function(e) data.table())
}
fs <- strsplit(opt$files, ",")[[1]]
keys <- rbindlist(lapply(fs, function(f) {
  d <- read_gz(f)
  if (!nrow(d) || !all(c("chr", "pos") %in% names(d))) return(data.table())
  unique(d[, .(chr = as.character(chr), pos = as.integer(pos))])
}), fill = TRUE)
# all-empty inputs -> rbindlist yields a 0-column data.table; write empty BED
if (!nrow(keys) || !all(c("chr", "pos") %in% names(keys))) { file.create(opt$bed); quit(status = 0) }
keys <- unique(keys[!is.na(chr) & !is.na(pos)])
if (!nrow(keys)) { file.create(opt$bed); quit(status = 0) }
setorder(keys, chr, pos)
fwrite(keys[, .(chr, start = pos - 1L, end = pos)], opt$bed, sep = "\t", col.names = FALSE)
