#!/usr/bin/env Rscript
# Convert a SNP list RData (data.frame `snplist`, column site = "chr:pos", the format of the
# original study's SNP_files) into the GRanges RData the pipeline's --snp_file expects:
# object snp_gr, width-1 ranges, unique and sorted.
#
# Usage: snp_list_to_granges.R <in.RData> <out.RData>
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })

a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 2)
e <- new.env()
load(a[1], envir = e)
x <- e[[ls(e)[1]]]
site <- if (is.data.frame(x)) x$site else as.character(x)

p <- tstrsplit(site, ":", fixed = TRUE)
snp_gr <- GRanges(seqnames = p[[1]], ranges = IRanges(start = as.integer(p[[2]]), width = 1))
snp_gr <- sort(unique(snp_gr))
save(snp_gr, file = a[2])
cat("input sites", length(site), "; unique positions", length(snp_gr), "; saved", a[2], "\n")
