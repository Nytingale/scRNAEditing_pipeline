#!/bin/bash
# ============================================================================
# Download mm39 reference genome and gene annotations for scRNAEditing pipeline
# ============================================================================
#
# Downloads:
#   1. mm39 FASTA from UCSC (GRCm39, UCSC chromosome naming)
#   2. NCBI RefSeq GTF from UCSC
#   3. Chrom sizes for extended BED generation
#
# Source: UCSC Genome Browser, matching original analysis (accessed 6 May 2024)
#
# Usage:
#   REFDIR=<dir> bash 01_download_genome_and_gtf.sh
#
# Output:
#   $REFDIR/mm39.fa          (indexed with samtools faidx + GATK dict)
#   $REFDIR/mm39.ncbiRefSeq.gtf
#   $REFDIR/mm39.chrom.sizes
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

echo "=== $(date) === Working directory: $REFDIR ==="

# ---- 1. Download mm39 FASTA ----

GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"
GENOME_GZ="mm39.fa.gz"
GENOME_FA="mm39.fa"

if [ -f "$GENOME_FA" ]; then
    echo "Genome already exists: $GENOME_FA"
else
    echo "=== $(date) === Downloading mm39 genome ==="
    wget --timestamping "$GENOME_URL" -O "$GENOME_GZ"
    echo "=== $(date) === Decompressing ==="
    zcat "$GENOME_GZ" > "$GENOME_FA"
    rm "$GENOME_GZ"
fi

if [ ! -f "${GENOME_FA}.fai" ]; then
    echo "=== $(date) === Indexing with samtools faidx ==="
    samtools faidx "$GENOME_FA"
fi

GENOME_DICT="${GENOME_FA%.fa}.dict"
if [ ! -f "$GENOME_DICT" ]; then
    echo "=== $(date) === Creating sequence dictionary ==="
    gatk CreateSequenceDictionary -R "$GENOME_FA"
fi

# ---- 2. Download NCBI RefSeq GTF ----

GTF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz"
GTF_GZ="mm39.ncbiRefSeq.gtf.gz"
GTF="mm39.ncbiRefSeq.gtf"

if [ -f "$GTF" ]; then
    echo "GTF already exists: $GTF"
else
    echo "=== $(date) === Downloading NCBI RefSeq GTF ==="
    wget --timestamping "$GTF_URL" -O "$GTF_GZ"
    echo "=== $(date) === Decompressing ==="
    gunzip "$GTF_GZ"
fi

# ---- 3. Download chrom sizes ----

CHROM_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes"
if [ ! -f "mm39.chrom.sizes" ]; then
    echo "=== $(date) === Downloading chrom sizes ==="
    wget "$CHROM_URL" -O "mm39.chrom.sizes"
fi

# ---- 4. Create main chromosome list ----

if [ ! -f "chrlist_main_chrs.txt" ]; then
    echo "=== $(date) === Creating main chromosome list ==="
    awk '$1 ~ /^chr[0-9]+$/ || $1 == "chrX" || $1 == "chrY"' mm39.chrom.sizes | \
        cut -f1 | sort -V > chrlist_main_chrs.txt
    echo "  Main chromosomes: $(wc -l < chrlist_main_chrs.txt)"
fi

# ---- 5. Verify ----

echo ""
echo "=== Verification ==="
echo "Genome: $(ls -lh "$GENOME_FA" | awk '{print $5}'), $(grep -c '^>' "$GENOME_FA") sequences"
echo "GTF: $(wc -l < "$GTF") lines"
echo "Chrom sizes: $(wc -l < mm39.chrom.sizes) sequences"
echo "Main chrs: $(wc -l < chrlist_main_chrs.txt)"

# ---- 6. Provenance ----

cat > download_provenance_genome.txt <<EOF
# Mouse reference genome and annotation provenance
# Generated: $(date -u +"%Y-%m-%d %H:%M:%S UTC")

Genome FASTA:
  File: $GENOME_FA
  Source: UCSC Genome Browser
  URL: $GENOME_URL
  Build: GRCm39 / mm39, UCSC chromosome naming (chr1, chr2, ...)
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
