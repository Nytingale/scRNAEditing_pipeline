#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
})

option_list = list(
  make_option("--input", type = "character",
              help = "Input SNP file: RData (data.frame with 'site' col as chr:pos) or BED (tab-delimited, cols 1=chr, 3=pos)"),
  make_option("--output", type = "character",
              help = "Output RData file containing a GRanges object named 'snp_gr'"),
  make_option("--format", type = "character", default = "auto",
              help = "Input format: 'rdata', 'bed', or 'auto' (detect by extension) [default: auto]")
)
opt = parse_args(OptionParser(option_list = option_list))

ext = tolower(tools::file_ext(opt$input))
fmt = opt$format
if (fmt == "auto") {
  fmt = if (ext %in% c("rdata", "rds")) "rdata" else "bed"
}
cat("Input:", opt$input, "\nDetected format:", fmt, "\n")

if (fmt == "rdata") {
  e = new.env()
  load(opt$input, envir = e)
  obj_names = ls(e)
  cat("Objects in RData:", paste(obj_names, collapse = ", "), "\n")

  snp_obj = NULL
  for (nm in obj_names) {
    obj = e[[nm]]
    if (inherits(obj, "GRanges")) {
      cat("Found GRanges object '", nm, "' with", length(obj), "ranges — using directly\n")
      snp_gr = obj
      save(snp_gr, file = opt$output)
      cat("Saved to", opt$output, "\n")
      quit(status = 0)
    }
    if (is.data.frame(obj) && "site" %in% colnames(obj)) {
      snp_obj = obj
      cat("Found data.frame '", nm, "' with", nrow(obj), "rows and 'site' column\n")
    }
  }

  if (is.null(snp_obj)) {
    stop("No data.frame with 'site' column or GRanges found in ", opt$input)
  }

  parts = data.table::tstrsplit(snp_obj$site, ":", fixed = TRUE)
  chr = parts[[1]]
  pos = as.integer(parts[[2]])

} else {
  snp_df = fread(opt$input, header = FALSE, select = c(1, 3))
  colnames(snp_df) = c("chr", "pos")
  chr = snp_df$chr
  pos = as.integer(snp_df$pos)
  cat("Read BED with", length(pos), "SNP positions\n")
}

valid = !is.na(pos)
cat("Valid positions:", sum(valid), "of", length(pos), "\n")
chr = chr[valid]
pos = pos[valid]

snp_gr = GRanges(seqnames = chr, ranges = IRanges(start = pos, end = pos))
snp_gr = sort(unique(snp_gr))
cat("GRanges object:", length(snp_gr), "unique SNP positions\n")

save(snp_gr, file = opt$output)
cat("Saved to", opt$output, "\n")
