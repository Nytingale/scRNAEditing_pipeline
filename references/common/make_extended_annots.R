#!/usr/bin/env Rscript
# Build the pipeline's --strand_annots / --alu_annots RData from the extended BEDs of
# generate_extended_bed.sh and (optionally) the Alu tables of human/05_alu_from_rmsk.sh.
# Objects: plus.annot, minus.annot (GRanges, mcols id + geneType) and, with --alu_*,
# alu.plus.annot, alu.minus.annot (mcols id). Strand is stored as "*": which strand an
# annotation belongs to is given by the object it is in.
#
# Usage: make_extended_annots.R --main_chrs chrlist_main_chrs.txt \
#          --plus <P>.plus.extended.sorted.bed --minus <P>.minus.extended.sorted.bed \
#          [--alu_plus alu_plus.txt --alu_minus alu_minus.txt] --out extended_annots_<genome>.RData
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1]
for (k in c("main_chrs", "plus", "minus", "out")) if (is.null(opt[[k]])) stop("missing --", k)
main_chrs <- readLines(opt$main_chrs)

# BED from gtf2bed: chr, start (0-based), end, id, score, strand, source, type, phase, attributes
bedToGRanges <- function(file) {
  df <- fread(file, nThread = 2, data.table = FALSE)[, -c(5, 7, 9, 10)]
  headers <- c("chr", "start", "end", "id", "strand", "geneType")
  colnames(df) <- headers[seq_len(length(df))]
  df <- df[df$chr %in% main_chrs, ]
  if ("strand" %in% colnames(df)) df$strand <- "*"
  df$start <- df$start + 1
  if (length(df) > 5) {
    with(df, GRanges(chr, IRanges(start, end), id = id, strand = strand, geneType = geneType))
  } else {
    with(df, GRanges(chr, IRanges(start, end), id = id, strand = strand))
  }
}

# Alu table: chr, start (0-based), end, repName, repFamily, strand
aluBedToGRanges <- function(file) {
  df <- fread(file, nThread = 2, data.table = FALSE,
              col.names = c("chr", "start", "end", "id", "drop", "strand"))
  df$drop <- NULL
  df <- df[df$chr %in% main_chrs, ]
  df$strand <- "*"
  df$start <- df$start + 1
  with(df, GRanges(chr, IRanges(start, end), id = id, strand = strand))
}

plus.annot <- bedToGRanges(opt$plus)
minus.annot <- bedToGRanges(opt$minus)
objs <- c("plus.annot", "minus.annot")
if (!is.null(opt$alu_plus)) {
  alu.plus.annot <- aluBedToGRanges(opt$alu_plus)
  alu.minus.annot <- aluBedToGRanges(opt$alu_minus)
  objs <- c(objs, "alu.plus.annot", "alu.minus.annot")
}
save(list = objs, file = opt$out)
for (o in objs) cat(o, length(get(o)), "\n")
cat("Saved", opt$out, "\n")
