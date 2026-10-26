#!/bin/bash
# ============================================================================
# Create SNP RData objects for the scRNAEditing pipeline (mouse)
# ============================================================================
#
# Creates two RData files from the extracted SNP coordinates:
#   1. GRCm39_snps_main_chrs.RData - string-based (legacy: snplist$site)
#   2. GRCm39_snps_GRanges.RData   - GRanges object (pipeline uses this)
#
# Prerequisite: Run 02_download_mgp_snp_vcf.sh first to create GRCm39_snp_coords.gz
#
# Usage:
#   REFDIR=<dir> bash 03_make_snp_rdata.sh
#   OR: Rscript 03_make_snp_rdata.R  (submitted via this wrapper)
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

echo "=== $(date) === Working directory: $REFDIR ==="

if [ ! -f "GRCm39_snp_coords.gz" ]; then
    echo "ERROR: GRCm39_snp_coords.gz not found. Run 02_download_mgp_snp_vcf.sh first."
    exit 1
fi

if [ ! -f "chrlist_main_chrs.txt" ]; then
    echo "ERROR: chrlist_main_chrs.txt not found. Run 01_download_genome_and_gtf.sh first."
    exit 1
fi

Rscript --vanilla - <<'RSCRIPT'
library(data.table)

cat("=== Loading SNP coordinates ===\n")
snplist <- fread(cmd = "zcat GRCm39_snp_coords.gz", nThread = 4, header = FALSE,
                 col.names = c("chr", "pos"))
cat(sprintf("  Total SNPs: %s\n", formatC(nrow(snplist), format = "d", big.mark = ",")))

main_chrs <- readLines("chrlist_main_chrs.txt")
cat(sprintf("  Main chromosomes: %d (%s ... %s)\n",
            length(main_chrs), main_chrs[1], tail(main_chrs, 1)))

# MGP VCF has been converted to UCSC names already (chr1, etc.)
snplist <- snplist[chr %in% main_chrs]
cat(sprintf("  SNPs on main chromosomes: %s\n",
            formatC(nrow(snplist), format = "d", big.mark = ",")))

# ---- String-based format (legacy) ----
cat("=== Creating string-based RData (legacy format) ===\n")
snplist_legacy <- data.frame(site = paste0(snplist$chr, ":", snplist$pos))
save(snplist_legacy, file = "GRCm39_snps_main_chrs.RData")
cat(sprintf("  Saved: GRCm39_snps_main_chrs.RData (%s)\n",
            format(file.size("GRCm39_snps_main_chrs.RData"), big.mark = ",")))

# ---- GRanges format (used by pipeline) ----
cat("=== Creating GRanges RData ===\n")
suppressPackageStartupMessages(library(GenomicRanges))

snp_gr <- GRanges(
    seqnames = snplist$chr,
    ranges = IRanges(start = snplist$pos, width = 1)
)
snp_gr <- sort(snp_gr)

cat(sprintf("  GRanges length: %s\n", formatC(length(snp_gr), format = "d", big.mark = ",")))
save(snp_gr, file = "GRCm39_snps_GRanges.RData")
cat(sprintf("  Saved: GRCm39_snps_GRanges.RData (%s)\n",
            format(file.size("GRCm39_snps_GRanges.RData"), big.mark = ",")))

cat("=== Done ===\n")
RSCRIPT

echo "=== $(date) === R processing complete ==="

# GRCm39_snp_coords.gz is kept: other analyses read it.
