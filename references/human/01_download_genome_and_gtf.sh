#!/bin/bash
# ============================================================================
# Download hg38 reference genome and gene annotations for scRNAEditing pipeline
# ============================================================================
#
# Downloads:
#   1. hg38 analysis set (no ALT contigs) from UCSC
#   2. NCBI RefSeq GTF from UCSC
#
# Reference genome:
#   Source: UCSC Genome Browser, hg38 analysis set (NO ALT contigs)
#   URL:    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/
#   File:   hg38.analysisSet.fa.gz  (NOT hg38.fullAnalysisSet which includes ALT scaffolds)
#   Build:  GRCh38, no ALT contigs, UCSC chromosome naming (chr1, chr2, ...)
#   Equivalent to NCBI: GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#   Note:   Non-canonical bases changed to N; repeats soft-masked (lower case)
#
# Gene annotation:
#   Source: UCSC Genome Browser, NCBI RefSeq track
#   URL:    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
#   Track:  ncbiRefSeq (RefSeq transcripts as aligned by NCBI)
#
# Usage:
#   REFDIR=<dir> bash 01_download_genome_and_gtf.sh
#
# Output:
#   $REFDIR/hg38.analysisSet.fa       (indexed with samtools faidx)
#   $REFDIR/hg38.ncbiRefSeq.gtf
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

echo "=== $(date) === Working directory: $REFDIR ==="

# ---- 1. Download hg38 analysis set FASTA ----

GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"
GENOME_GZ="hg38.analysisSet.fa.gz"
GENOME_FA="hg38.analysisSet.fa"

if [ -f "$GENOME_FA" ]; then
    echo "Genome already exists: $GENOME_FA"
else
    echo "=== $(date) === Downloading hg38 analysis set ==="
    wget --timestamping "$GENOME_URL" -O "$GENOME_GZ"

    echo "=== $(date) === Decompressing ==="
    gunzip "$GENOME_GZ"
fi

# Index with samtools
if [ ! -f "${GENOME_FA}.fai" ]; then
    echo "=== $(date) === Indexing with samtools faidx ==="
    samtools faidx "$GENOME_FA"
fi

# Create GATK sequence dictionary
GENOME_DICT="${GENOME_FA%.fa}.dict"
if [ ! -f "$GENOME_DICT" ]; then
    echo "=== $(date) === Creating sequence dictionary ==="
    gatk CreateSequenceDictionary -R "$GENOME_FA"
fi

# ---- 2. Download NCBI RefSeq GTF ----

GTF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
GTF_GZ="hg38.ncbiRefSeq.gtf.gz"
GTF="hg38.ncbiRefSeq.gtf"

if [ -f "$GTF" ]; then
    echo "GTF already exists: $GTF"
else
    echo "=== $(date) === Downloading NCBI RefSeq GTF ==="
    wget --timestamping "$GTF_URL" -O "$GTF_GZ"

    echo "=== $(date) === Decompressing ==="
    gunzip "$GTF_GZ"
fi

# ---- 3. Chromosome sizes and main chromosome list (chr1-22, X, Y) ----

cut -f1,2 "${GENOME_FA}.fai" > hg38.analysisSet.chrom.sizes
if [ ! -f "chrlist_main_chrs.txt" ]; then
    awk '$1 ~ /^chr[0-9]+$/ || $1 == "chrX" || $1 == "chrY"' hg38.analysisSet.chrom.sizes | \
        cut -f1 | sort -V > chrlist_main_chrs.txt
fi

# ---- 4. Verify downloads ----

echo ""
echo "=== Verification ==="
echo "Genome FASTA:"
head -1 "$GENOME_FA"
echo "  Chromosomes: $(grep -c '^>' "$GENOME_FA")"
echo "  Size: $(ls -lh "$GENOME_FA" | awk '{print $5}')"
echo "  Index: $(wc -l < "${GENOME_FA}.fai") sequences in .fai"

echo ""
echo "RefSeq GTF:"
echo "  Lines: $(wc -l < "$GTF")"
echo "  Chromosomes: $(awk -F'\t' '$1 !~ /^#/' "$GTF" | cut -f1 | sort -u | head -5) ..."
echo "  Size: $(ls -lh "$GTF" | awk '{print $5}')"

# ---- 5. Record provenance ----

cat > download_provenance_genome.txt <<EOF
# Reference genome and annotation provenance
# Generated: $(date -u +"%Y-%m-%d %H:%M:%S UTC")

Genome FASTA:
  File: $GENOME_FA
  Source: UCSC Genome Browser
  URL: $GENOME_URL
  Build: GRCh38 / hg38 analysis set (no ALT contigs)
  Description: Primary assembly with non-canonical bases changed to N
  md5: $(md5sum "$GENOME_FA" | awk '{print $1}')

Gene Annotation:
  File: $GTF
  Source: UCSC Genome Browser, NCBI RefSeq track
  URL: $GTF_URL
  Track: ncbiRefSeq
  md5: $(md5sum "$GTF" | awk '{print $1}')
EOF

echo ""
echo "=== $(date) === Done ==="
echo "Provenance recorded in: download_provenance_genome.txt"
