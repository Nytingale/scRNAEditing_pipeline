#!/bin/bash
# ============================================================================
# STARsolo alignment for teratoma (PM) dataset
# ============================================================================
# 8 SRR accessions → align per-SRR → merge replicates per teratoma
# (GEO GSE248941 / GSM7923327–330; pairing corrected 2026-09-05):
#   T1: SRR26992458 (T1_1) + SRR26992459 (T1_2)
#   T2: SRR26992455 (T2_1) + SRR26992456 (T2_2)
#   T3: SRR26992452 (T3_1) + SRR26992453 (T3_2)
#   T4: SRR26992449 (T4_1) + SRR26992450 (T4_2)
#
# Platform: 10x v3 (CB=16bp, UMI=12bp, whitelist=3M-february-2018.txt)
# Read length: 101bp → sjdbOverhang=100
# R1 = barcode read, R2 = cDNA read
# STARsolo expects: --readFilesIn <cDNA> <barcode> (R2 R1)
# ============================================================================

set -euo pipefail

SAMTOOLS="${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/samtools"
if [ ! -x "$SAMTOOLS" ]; then
    echo "ERROR: samtools not found in scRNAEditing env: $SAMTOOLS" >&2
    exit 1
fi

REFDIR="${SCRNA_DIR:-.}/human/ref_annots"
FQDIR="${SCRNA_DIR:-.}/human/pm/fastq_files"
OUTDIR="${SCRNA_DIR:-.}/human/pm/01_starsolobams"
STARDB="${REFDIR}/stardb_hg38_100"
WHITELIST="${SCRNA_DIR:-.}/nextflow/assets/whitelists/3M-february-2018.txt"
NTHREADS=${NTHREADS:-8}

SAMPLE=${1:?Usage: bash align_teratoma.sh <SRR_accession>}

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== $(date) === Aligning $SAMPLE ==="

R1="${FQDIR}/${SAMPLE}_1.fastq.gz"
R2="${FQDIR}/${SAMPLE}_2.fastq.gz"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: FASTQ files not found for $SAMPLE"
    exit 1
fi

STAR \
    --genomeDir "$STARDB" \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBwhitelist "$WHITELIST" \
    --readFilesIn "$R2" "$R1" \
    --readFilesCommand zcat \
    --soloFeatures GeneFull SJ \
    --twopassMode Basic \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN MD \
    --outFilterType BySJout \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM 40000000000 \
    --genomeLoad NoSharedMemory \
    --outFileNamePrefix "${SAMPLE}_" \
    --runThreadN "$NTHREADS"

echo "=== $(date) === Indexing BAM ==="
echo "Using $($SAMTOOLS --version | head -1) [$SAMTOOLS]"
"$SAMTOOLS" index "${SAMPLE}_Aligned.sortedByCoord.out.bam"

echo "=== $(date) === Done: $SAMPLE ==="
