#!/bin/bash
# Extended gene annotation BEDs for stranding sites (exons, CDS, UTRs, introns, intergenic),
# one file per strand. Run in a directory containing:
#   <PREFIX>.<ANNOTTYPE>.gtf.gz   gene annotation
#   <PREFIX>.chrom.sizes          chromosome sizes (e.g. cut -f1,2 genome.fa.fai)
# Writes <PREFIX>.<ANNOTTYPE>.{plus,minus}.extended.sorted.bed and <PREFIX>_geneTypes.txt.
# Requires gtf2bed (BEDOPS) and bedtools.
#
# Usage: generate_extended_bed.sh <PREFIX> <ANNOTTYPE>     e.g. GRCh38_analysis_set ncbiRefSeq
set -eu

PREFIX=$1
ANNOTTYPE=$2
P="${PREFIX}.${ANNOTTYPE}"

# Main chromosomes only: drop fix / alt / random / unplaced contigs and the mitochondrion.
# convert2bed (BEDOPS 2.4.39-2.4.42) can abort with heap corruption after writing all records,
# so its exit status is not used: the record count is checked instead.
main_gtf() {
    zcat "${P}.gtf.gz" | \
    grep -v -E "^chr[a-zA-Z0-9]+_[a-zA-Z0-9]+_fix" | \
    grep -v -E "^chr[a-zA-Z0-9]+_[a-zA-Z0-9]+_alt" | \
    grep -v -E "^chr[a-zA-Z0-9]+_[a-zA-Z0-9]+_random" | \
    grep -v "^chrUn" | grep -v "^chrM" | grep -v "^MT" | grep -v "^Sca"
}
( ulimit -c 0; main_gtf | gtf2bed --do-not-sort 2>/dev/null | sort-bed --max-mem 500M - | gzip -f > "${P}.bed.gz" )
n_gtf=$(main_gtf | grep -vc '^#' || true)
n_bed=$(zcat "${P}.bed.gz" | wc -l)
[ "$n_gtf" -eq "$n_bed" ] || { echo "gtf2bed wrote $n_bed of $n_gtf records" >&2; exit 1; }

zcat "${P}.bed.gz" | cut -f8 | sort | uniq -c | sort -k1rn > "${PREFIX}_geneTypes.txt"

sort -k1,1 -k2,2n "${PREFIX}.chrom.sizes" > "${PREFIX}.genome"

zcat "${P}.bed.gz" | awk 'BEGIN{OFS="\t";} $6=="+" {print $0}' | gzip -f > "${P}.plus.bed.gz"
zcat "${P}.bed.gz" | awk 'BEGIN{OFS="\t";} $6=="-" {print $0}' | gzip -f > "${P}.minus.bed.gz"

# Source name written in the intergenic rows (GTF column 2, e.g. ncbiRefSeq.2022-10-28)
annotType=$(zcat "${P}.gtf.gz" | grep -v "#" | cut -f2 | head -n1 || true)

generate_extended_bed(){
    STRAND=$1

    zcat "${P}.${STRAND}.bed.gz" | \
    awk 'BEGIN{FS="\t"; OFS="\t";} $8=="exon" {print $0}' | \
    bedtools sort | bedtools merge -i - | gzip -f > "${P}.${STRAND}_exon_merged.bed.gz"

    # Introns: transcript / gene bodies minus merged exons
    zcat "${P}.${STRAND}.bed.gz" | \
    awk 'BEGIN{FS="\t"; OFS="\t";} $8=="transcript" || $8=="gene" {print $0}' | \
    bedtools sort | bedtools subtract -a stdin -b "${P}.${STRAND}_exon_merged.bed.gz" | \
    awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,"intron",$9,$10}' | \
    gzip -f > "${P}.${STRAND}_intron.bed.gz"

    # Intergenic: genome not covered by a transcript on this strand
    zcat "${P}.${STRAND}.bed.gz" | \
    awk 'BEGIN{OFS="\t";} $8=="transcript" || $8=="gene" {print $0}' | \
    bedtools sort | bedtools complement -i stdin -g "${PREFIX}.genome" | \
    awk -v annotType="$annotType" 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,"NA","NA","NA",annotType,"intergenic","NA","NA"}' | \
    gzip -f > "${P}.${STRAND}_intergenic.bed.gz"

    cat "${P}.${STRAND}.bed.gz" "${P}.${STRAND}_intron.bed.gz" "${P}.${STRAND}_intergenic.bed.gz" \
        > "${P}.${STRAND}.extended.bed.gz"
    zcat "${P}.${STRAND}.extended.bed.gz" | bedtools sort > "${P}.${STRAND}.extended.sorted.bed"

    rm -f "${P}.${STRAND}_exon_merged.bed.gz" "${P}.${STRAND}_intron.bed.gz" \
          "${P}.${STRAND}_intergenic.bed.gz" "${P}.${STRAND}.extended.bed.gz" "${P}.${STRAND}.bed.gz"
}

echo "Processing plus strand..."
generate_extended_bed "plus"
echo "Processing minus strand..."
generate_extended_bed "minus"

rm -f "${PREFIX}.genome" "${P}.bed.gz"
