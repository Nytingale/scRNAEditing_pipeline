#!/bin/bash
# ============================================================================
# Build STAR genome index for mm39 (mouse)
# ============================================================================
#
# Builds a STAR index with sjdbOverhang=49 (matching the 50bp cDNA reads
# from JZL 10x v2 mouse cortex datasets).
#
# If additional mouse datasets with different read lengths are processed,
# additional indices can be built by passing a different overhang value.
#
# Prerequisites:
#   - Run 01_download_genome_and_gtf.sh first
#   - STAR must be available in PATH
#
# Usage:
#   REFDIR=<dir> bash 04_build_star_index.sh [OVERHANG]
#   Default OVERHANG=49
#
# Output:
#   $REFDIR/stardb_mm39_<overhang>/
#
# Resource requirements:
#   ~32 GB RAM, ~12 CPU cores, ~1-2 hours
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

GENOME_FA="mm39.fa"
GTF="mm39.ncbiRefSeq.gtf"
NTHREADS=${NTHREADS:-12}
OVERHANG=${1:-49}

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

OUTDIR="stardb_mm39_${OVERHANG}"

if [ -f "${OUTDIR}/SA" ]; then
    echo "STAR index already exists: $OUTDIR"
    exit 0
fi

echo "=== $(date) === Building STAR index (sjdbOverhang=${OVERHANG}) ==="
mkdir -p "$OUTDIR"

STAR \
    --runMode genomeGenerate \
    --genomeDir "$OUTDIR" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang "$OVERHANG" \
    --runThreadN "$NTHREADS"

echo "=== $(date) === Index built: $OUTDIR ==="
echo "  Size: $(du -sh "$OUTDIR" | cut -f1)"
