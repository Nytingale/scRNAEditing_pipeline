#!/bin/bash
# extended_annots_mm39.RData (--strand_annots): extended RefSeq gene annotation per strand.
# No Alu annotation for mouse (Alu elements are primate-specific). Needs 01 (GTF, chrom
# sizes, chromosome list). Requires gtf2bed (BEDOPS), bedtools and R with GenomicRanges, data.table.
#
# Usage: REFDIR=<dir> bash 05_extended_annots.sh
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REFDIR"

P=mm39
mkdir -p extended_bed && cd extended_bed
gzip -c ../mm39.ncbiRefSeq.gtf > ${P}.ncbiRefSeq.gtf.gz
cp ../mm39.chrom.sizes ${P}.chrom.sizes
bash "$HERE/../common/generate_extended_bed.sh" "$P" ncbiRefSeq
cd ..

Rscript "$HERE/../common/make_extended_annots.R" \
    --main_chrs chrlist_main_chrs.txt \
    --plus  extended_bed/${P}.ncbiRefSeq.plus.extended.sorted.bed \
    --minus extended_bed/${P}.ncbiRefSeq.minus.extended.sorted.bed \
    --out extended_annots_${P}.RData
