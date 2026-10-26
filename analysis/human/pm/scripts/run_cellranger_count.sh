#!/bin/bash
# CellRanger count on ADAR-KO teratomas (T1-T4)
# SGE array job: -t 1-4
# Depends on: build_cellranger_ref.sh completing first
# Date: 2026-09-02

set -euo pipefail

CELLRANGER="${SOFTWARE_DIR:-/opt/software}/cellranger-10.0.0/bin/cellranger"
CRREF="${SCRNA_DIR:-.}/human/pm/cellranger_ref/GRCh38_and_GRCm38"
CRDIR="${SCRNA_DIR:-.}/human/pm/cellranger_fqs"
OUTDIR="${SCRNA_DIR:-.}/human/pm/cellranger_out"

SAMPLES=("T1" "T2" "T3" "T4")
SAMPLE="${SAMPLES[$SGE_TASK_ID-1]}"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Remove previous run directory if it exists (CellRanger won't overwrite)
rm -rf "$OUTDIR/$SAMPLE"

echo "=== Running cellranger count for $SAMPLE ==="
$CELLRANGER count \
  --id="$SAMPLE" \
  --transcriptome="$CRREF" \
  --fastqs="$CRDIR" \
  --sample="$SAMPLE" \
  --localcores=16 \
  --localmem=60 \
  --create-bam=false

echo "=== $SAMPLE complete ==="
echo "Metrics:"
cat "$OUTDIR/$SAMPLE/outs/metrics_summary.csv"
