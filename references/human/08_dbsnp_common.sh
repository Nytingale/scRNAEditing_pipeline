#!/bin/bash
# Common SNPs of dbSNP build 156 for --snp_file: every record of the dbSNP 156 VCF (from
# 02_download_dbsnp_vcf.sh) with INFO/COMMON=1, all variant classes (SNV, indel, MNV),
# one entry per VCF record. The intermediate dbSNP156_common.RData (data.frame snplist,
# column site = "chr:pos") is the original study's list, identical in content and order.
#
# Usage: REFDIR=<dir> bash 08_dbsnp_common.sh
# Needs bcftools and R (GenomicRanges, data.table); about 20 min and 32 GB memory.
# Output: $REFDIR/dbSNP156_common.RData, $REFDIR/dbSNP156_common_GRanges.RData
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REFDIR"
[ -f dbSNP156_GRCh38.vcf.gz ] || { echo "run 02_download_dbsnp_vcf.sh first" >&2; exit 1; }

bcftools query -i 'INFO/COMMON=1' -f '%CHROM\t%POS\n' dbSNP156_GRCh38.vcf.gz \
    | gzip -1 > dbSNP156_common.tsv.gz

Rscript -e '
suppressPackageStartupMessages(library(data.table))
x <- fread(cmd = "zcat dbSNP156_common.tsv.gz", header = FALSE, col.names = c("chr", "pos"))
snplist <- data.frame(site = paste0(x$chr, ":", x$pos))
save(snplist, file = "dbSNP156_common.RData")
cat(nrow(snplist), "common records\n")'

Rscript "$HERE/../common/snp_list_to_granges.R" dbSNP156_common.RData dbSNP156_common_GRanges.RData
rm -f dbSNP156_common.tsv.gz
