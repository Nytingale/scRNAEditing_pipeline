#!/bin/bash
# extended_annots_GRCh38_analysis_set.RData (--strand_annots and --alu_annots): extended
# RefSeq gene annotation per strand plus Alu elements. Needs 01 (FASTA index, chromosome list)
# and 05 (Alu tables). Requires gtf2bed (BEDOPS), bedtools and R with GenomicRanges, data.table.
#
# The gene annotation is NCBI's RefSeq GTF for the GRCh38 analysis set (UCSC chromosome
# names, with gene rows), not the UCSC GTF given to --genome_gtf. NCBI replaces this file
# with each annotation release; the file used for the paper (2024) is on Zenodo.
#
# Usage: REFDIR=<dir> bash 06_extended_annots.sh
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REFDIR"

P=GRCh38_analysis_set
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
mkdir -p extended_bed && cd extended_bed
[ -f ${P}.ncbiRefSeq.gtf.gz ] || wget -O ${P}.ncbiRefSeq.gtf.gz "$URL"
cut -f1,2 ../hg38.analysisSet.fa.fai > ${P}.chrom.sizes
bash "$HERE/../common/generate_extended_bed.sh" "$P" ncbiRefSeq
cd ..

Rscript "$HERE/../common/make_extended_annots.R" \
    --main_chrs chrlist_main_chrs.txt \
    --plus  extended_bed/${P}.ncbiRefSeq.plus.extended.sorted.bed \
    --minus extended_bed/${P}.ncbiRefSeq.minus.extended.sorted.bed \
    --alu_plus hg38_alu_plus.txt --alu_minus hg38_alu_minus.txt \
    --out extended_annots_${P}.RData
