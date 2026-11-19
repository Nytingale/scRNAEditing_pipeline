#!/bin/bash
# ============================================================================
# Build the STAR index for Octopus bimaculoides (obim.fa + obim.ncbiRefSeq.gtf from 01)
# ============================================================================
#
# sjdbOverhang = cDNA read length - 1: 132 for the optic-lobe 10x data (PRJNA854179, 133-nt cDNA reads), as in
# the original study. The genome has ~145,000 scaffolds, so genomeGenerate needs ~120 GB RAM.
#
# Usage:
#   REFDIR=<dir> bash 02_build_star_index.sh [OVERHANG]      (default 132)
#
# Output:
#   $REFDIR/stardb_obim_<overhang>/
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
cd "$REFDIR"
OVERHANG=${1:-132}
NTHREADS=${NTHREADS:-12}
OUT="stardb_obim_${OVERHANG}"

for f in obim.fa obim.ncbiRefSeq.gtf; do
    [ -f "$f" ] || { echo "ERROR: $REFDIR/$f not found; run 01_download_genome_and_gtf.sh first" >&2; exit 1; }
done
if [ -f "${OUT}/SA" ]; then echo "STAR index already exists: $OUT"; exit 0; fi

echo "=== $(date) === Building STAR index (sjdbOverhang=${OVERHANG}) ==="
mkdir -p "$OUT"
STAR \
    --runMode genomeGenerate \
    --genomeDir "$OUT" \
    --genomeFastaFiles obim.fa \
    --sjdbGTFfile obim.ncbiRefSeq.gtf \
    --sjdbOverhang "$OVERHANG" \
    --limitGenomeGenerateRAM 120000000000 \
    --outTmpDir "${REFDIR}/_STARtmp_${OVERHANG}" \
    --runThreadN "$NTHREADS"
echo "=== $(date) === Index built: $OUT ($(du -sh "$OUT" | cut -f1)) ==="
