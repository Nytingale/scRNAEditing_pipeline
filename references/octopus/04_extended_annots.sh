#!/bin/bash
# extended_annots_obim.RData (--strand_annots): extended RefSeq gene annotation per strand, chr1-chr30.
# No Alu annotation (Alu elements are primate-specific). Needs 01 (GTF, chrom sizes, chromosome list).
# Requires gtf2bed (BEDOPS), bedtools and R with GenomicRanges, data.table (the refbuild environment).
# The original study's file is annotations/annotation_files/extended_annots_obim.RData; this rebuilds it.
#
# Usage: REFDIR=<dir> bash 04_extended_annots.sh
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REFDIR"

P=obim
mkdir -p extended_bed && cd extended_bed
gzip -c ../obim.ncbiRefSeq.gtf > ${P}.ncbiRefSeq.gtf.gz
cp ../obim.chrom.sizes ${P}.chrom.sizes
bash "$HERE/../common/generate_extended_bed.sh" "$P" ncbiRefSeq
cd ..

Rscript "$HERE/../common/make_extended_annots.R" \
    --main_chrs chrlist_main_chrs.txt \
    --plus  extended_bed/${P}.ncbiRefSeq.plus.extended.sorted.bed \
    --minus extended_bed/${P}.ncbiRefSeq.minus.extended.sorted.bed \
    --out extended_annots_${P}.RData
