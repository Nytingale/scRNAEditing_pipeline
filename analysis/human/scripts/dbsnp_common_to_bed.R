#!/usr/bin/env Rscript
# =============================================================================
# dbsnp_common_to_bed.R
#
# Export the dbSNP156 common SNP GRanges (object `snp_gr`, the same file and
# object nextflow/bin/downstream_filter.R uses for snp_subtract=on) to a BED3
# of single positions, for bedtools-based SNP flagging in GT v6.
#
# The pipeline matches SNPs with findOverlaps(type = "equal") against width-1
# sites, so only width-1 ranges can ever match; wider ranges are dropped here
# (count reported) to keep that behaviour identical.
#
# Input : annotations/SNP_files/dbSNP156_common_GRanges.RData
# Output: human/pm/gt_eval_v6/snp/dbSNP156_common.bed   (chr, pos-1, pos)
#
# Usage : Rscript dbsnp_common_to_bed.R <in.RData> <out.bed>
# =============================================================================

suppressPackageStartupMessages(library(GenomicRanges))

a = commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 2)

e = new.env()
load(a[1], envir = e)
if (!"snp_gr" %in% ls(e)) stop("object 'snp_gr' not found in ", a[1])
gr = e$snp_gr
stopifnot(inherits(gr, "GRanges"))

chr = as.character(seqnames(gr))
chr[chr %in% c("MT", "M")] = "chrM"
chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
w1 = width(gr) == 1

bed = unique(data.frame(chr = chr[w1], start = start(gr)[w1] - 1L, end = start(gr)[w1]))
bed = bed[order(bed$chr, bed$start), ]
write.table(bed, a[2], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("snp_gr ranges       :", length(gr), "\n")
cat("width-1 kept        :", sum(w1), "\n")
cat("wider ranges dropped:", sum(!w1), "\n")
cat("unique BED positions:", nrow(bed), "\n")
cat("chromosomes         :", paste(head(unique(bed$chr), 30), collapse = " "), "\n")
