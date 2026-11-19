#!/bin/bash
# ============================================================================
# Download the Octopus bimaculoides reference (NCBI RefSeq GCF_001194135.2, ASM119413v2) for scRNAEditing
# ============================================================================
#
# Same source and renaming as the original study (accessed 21 Jun 2024): the 30 chromosomes and the mitochondrion
# are renamed from RefSeq accessions to chr<n> / chrMT with the NCBI assembly report; the unplaced NW_ scaffolds
# keep their accessions and stay in the FASTA (so reads from them are not forced onto the chromosomes), but are
# left out of obim.chrom.sizes and the main-chromosome list, as are chrMT. These names match
# extended_annots_obim.RData and obim_SNPs_GRanges.RData (chr1-chr30).
# FASTA headers are cut to the sequence name.
#
# Usage:
#   REFDIR=<dir> bash 01_download_genome_and_gtf.sh
#
# Output:
#   $REFDIR/obim.fa (+ .fai, obim.dict)
#   $REFDIR/obim.ncbiRefSeq.gtf          chromosomes renamed
#   $REFDIR/obim.chrom.sizes             chr1-chr30
#   $REFDIR/chrlist_main_chrs.txt        chr1-chr30
#   $REFDIR/chr_name_conv_obim.txt       RefSeq accession -> chr name
#   $REFDIR/GCF_001194135.2_ASM119413v2_assembly_report.txt
#   $REFDIR/download_provenance_genome.txt
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

ACC=GCF_001194135.2_ASM119413v2
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/194/135/${ACC}"
echo "=== $(date) === Working directory: $REFDIR ==="

wget -q -O md5checksums.ncbi.txt "${URL}/md5checksums.txt"
for f in "${ACC}_genomic.fna.gz" "${ACC}_genomic.gtf.gz" "${ACC}_assembly_report.txt"; do
    if [ ! -f "$f" ]; then
        echo "=== $(date) === Downloading $f ==="
        wget -q -O "$f" "${URL}/$f"
    fi
    want=$(awk -v f="./$f" '$2 == f {print $1}' md5checksums.ncbi.txt)
    got=$(md5sum "$f" | cut -d' ' -f1)
    [ "$want" = "$got" ] || { echo "ERROR: md5 mismatch for $f ($got vs NCBI $want)" >&2; exit 1; }
done

# Assembled molecules: RefSeq accession (col 7) -> chr<Assigned-Molecule> (col 3); MT -> chrMT
grep -v '^#' "${ACC}_assembly_report.txt" | \
    awk -F'\t' '$2 == "assembled-molecule" {print $7 "\tchr" $3}' > chr_name_conv_obim.txt
[ "$(wc -l < chr_name_conv_obim.txt)" -eq 31 ] || { echo "ERROR: expected 30 chromosomes + MT" >&2; exit 1; }

if [ ! -f obim.fa ]; then
    echo "=== $(date) === Renaming FASTA ==="
    zcat "${ACC}_genomic.fna.gz" | awk -F'\t' 'NR == FNR {m[$1] = $2; next}
        /^>/ {n = substr($1, 2); sub(/ .*/, "", n); print ">" (n in m ? m[n] : n); next} {print}' \
        chr_name_conv_obim.txt - > obim.fa.tmp
    mv obim.fa.tmp obim.fa
fi
[ -f obim.fa.fai ] || samtools faidx obim.fa
[ -f obim.dict ] || gatk CreateSequenceDictionary -R obim.fa -O obim.dict

if [ ! -f obim.ncbiRefSeq.gtf ]; then
    echo "=== $(date) === Renaming GTF ==="
    zcat "${ACC}_genomic.gtf.gz" | awk 'BEGIN {FS = OFS = "\t"} NR == FNR {m[$1] = $2; next}
        /^#/ {print; next} {if ($1 in m) $1 = m[$1]; print}' chr_name_conv_obim.txt - > obim.ncbiRefSeq.gtf.tmp
    mv obim.ncbiRefSeq.gtf.tmp obim.ncbiRefSeq.gtf
fi

grep -E '^chr[0-9]+\b' obim.fa.fai | cut -f1,2 | sort -k1,1V > obim.chrom.sizes
cut -f1 obim.chrom.sizes > chrlist_main_chrs.txt
[ "$(wc -l < chrlist_main_chrs.txt)" -eq 30 ] || { echo "ERROR: expected chr1-chr30" >&2; exit 1; }

{
    echo "# Reference genome and annotation provenance"
    echo "# Generated: $(date -u '+%F %T') UTC"
    echo
    echo "Genome FASTA:"
    echo "  File: obim.fa (chromosomes renamed with chr_name_conv_obim.txt; headers cut to the name)"
    echo "  Source: NCBI RefSeq ${ACC}"
    echo "  URL: ${URL}/${ACC}_genomic.fna.gz"
    echo "  md5 (download): $(md5sum "${ACC}_genomic.fna.gz" | cut -d' ' -f1)"
    echo "  md5 (obim.fa): $(md5sum obim.fa | cut -d' ' -f1)"
    echo "  Sequences: $(wc -l < obim.fa.fai) (30 chromosomes, chrMT, $(grep -c '^NW_' obim.fa.fai) NW_ scaffolds)"
    echo
    echo "Gene Annotation:"
    echo "  File: obim.ncbiRefSeq.gtf (chromosomes renamed)"
    echo "  URL: ${URL}/${ACC}_genomic.gtf.gz"
    echo "  md5 (download): $(md5sum "${ACC}_genomic.gtf.gz" | cut -d' ' -f1)"
    echo "  md5 (obim.ncbiRefSeq.gtf): $(md5sum obim.ncbiRefSeq.gtf | cut -d' ' -f1)"
    echo
    echo "Assembly report: ${URL}/${ACC}_assembly_report.txt"
    echo "  md5: $(md5sum "${ACC}_assembly_report.txt" | cut -d' ' -f1)"
} > download_provenance_genome.txt
rm -f "${ACC}_genomic.fna.gz"
echo "=== $(date) === Done ==="
