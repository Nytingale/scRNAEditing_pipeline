#!/bin/bash
# ============================================================================
# Build STAR genome indices for scRNAEditing pipeline
# ============================================================================
#
# Builds two STAR indices with different sjdbOverhang values:
#   - sjdbOverhang=100: For datasets with longer cDNA reads (e.g., teratoma)
#   - sjdbOverhang=49:  For 10x v3 datasets (R2 cDNA read is ~50bp, so 50-1=49)
#
# sjdbOverhang should be set to (read length - 1) for optimal splice junction
# detection. Different datasets may have different read lengths.
#
# Prerequisites:
#   - Run 01_download_genome_and_gtf.sh first
#   - STAR must be available in PATH
#
# Usage:
#   REFDIR=<dir> bash 03_build_star_index.sh [OVERHANG]
#
#   If OVERHANG is provided, only that index is built.
#   If omitted, both 100 and 49 are built.
#
# Output:
#   $REFDIR/stardb_hg38_100/   (sjdbOverhang=100)
#   $REFDIR/stardb_hg38_49/    (sjdbOverhang=49)
#
# Resource requirements:
#   ~32 GB RAM, ~12 CPU cores, ~1-2 hours per index
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

GENOME_FA="hg38.analysisSet.fa"
GTF="hg38.ncbiRefSeq.gtf"
NTHREADS=${NTHREADS:-12}

if [ ! -f "$GENOME_FA" ]; then
    echo "ERROR: Genome FASTA not found: $REFDIR/$GENOME_FA"
    echo "Run 01_download_genome_and_gtf.sh first."
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF not found: $REFDIR/$GTF"
    echo "Run 01_download_genome_and_gtf.sh first."
    exit 1
fi

build_index() {
    local overhang=$1
    local outdir="stardb_hg38_${overhang}"

    if [ -f "${outdir}/SA" ]; then
        echo "STAR index already exists: $outdir"
        return 0
    fi

    echo "=== $(date) === Building STAR index (sjdbOverhang=${overhang}) ==="
    mkdir -p "$outdir"

    STAR \
        --runMode genomeGenerate \
        --genomeDir "$outdir" \
        --genomeFastaFiles "$GENOME_FA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang "$overhang" \
        --runThreadN "$NTHREADS"

    echo "=== $(date) === Index built: $outdir ==="
    echo "  Size: $(du -sh "$outdir" | cut -f1)"
}

if [ $# -ge 1 ]; then
    build_index "$1"
else
    build_index 100
    build_index 49
fi

echo "=== $(date) === All STAR indices complete ==="
