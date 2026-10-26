#!/bin/bash
# CellRanger mkgtf + mkref: build dual hg38/mm10 reference
# Depends on: download_cellranger_refs.sh completing first
# Date: 2026-09-02

set -euo pipefail

REFDIR="${SCRNA_DIR:-.}/human/pm/cellranger_ref"
CELLRANGER="${SOFTWARE_DIR:-/opt/software}/cellranger-10.0.0/bin/cellranger"

cd "$REFDIR"

echo "=== Filtering GTFs (protein_coding only) ==="
$CELLRANGER mkgtf \
  "$REFDIR/Homo_sapiens.GRCh38.98.gtf" \
  "$REFDIR/Homo_sapiens.GRCh38.98.filtered.gtf" \
  --attribute=gene_biotype:protein_coding

$CELLRANGER mkgtf \
  "$REFDIR/Mus_musculus.GRCm38.98.gtf" \
  "$REFDIR/Mus_musculus.GRCm38.98.filtered.gtf" \
  --attribute=gene_biotype:protein_coding

echo "=== Building dual-species reference ==="
$CELLRANGER mkref \
  --genome=GRCh38 \
  --fasta="$REFDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --genes="$REFDIR/Homo_sapiens.GRCh38.98.filtered.gtf" \
  --genome=GRCm38 \
  --fasta="$REFDIR/Mus_musculus.GRCm38.dna.primary_assembly.fa" \
  --genes="$REFDIR/Mus_musculus.GRCm38.98.filtered.gtf" \
  --nthreads=16 \
  --memgb=100

echo "=== mkref complete ==="
ls -lh "$REFDIR/GRCh38_and_GRCm38/"
