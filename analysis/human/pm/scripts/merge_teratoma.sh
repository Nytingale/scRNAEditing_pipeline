#!/bin/bash
# ============================================================================
# Merge teratoma BAMs: per-SRR BAMs → per-teratoma merged BAM
# ============================================================================
# GEO GSE248941 / GSM7923327–330 (confirmed 2026-09-05):
#   T1: SRR26992458 (T1_1) + SRR26992459 (T1_2)  GSM7923327
#   T2: SRR26992455 (T2_1) + SRR26992456 (T2_2)  GSM7923328
#   T3: SRR26992452 (T3_1) + SRR26992453 (T3_2)  GSM7923329
#   T4: SRR26992449 (T4_1) + SRR26992450 (T4_2)  GSM7923330
# Each GSM is one large 10x lane + one small second run of the same library.
# ============================================================================

set -euo pipefail

SAMTOOLS="${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/samtools"
if [ ! -x "$SAMTOOLS" ]; then
    echo "ERROR: samtools not found in scRNAEditing env: $SAMTOOLS" >&2
    exit 1
fi
echo "Using $($SAMTOOLS --version | head -1) [$SAMTOOLS]"

BAMDIR="${SCRNA_DIR:-.}/human/pm/01_starsolobams"
cd "$BAMDIR"

declare -A TERATOMA_SRRS
TERATOMA_SRRS[T1]="SRR26992458 SRR26992459"
TERATOMA_SRRS[T2]="SRR26992455 SRR26992456"
TERATOMA_SRRS[T3]="SRR26992452 SRR26992453"
TERATOMA_SRRS[T4]="SRR26992449 SRR26992450"

for TERA in T1 T2 T3 T4; do
    SRRS=(${TERATOMA_SRRS[$TERA]})
    INPUTS=""
    ALL_EXIST=true

    for SRR in "${SRRS[@]}"; do
        BAM="${SRR}_Aligned.sortedByCoord.out.bam"
        if [ ! -f "$BAM" ]; then
            echo "WARNING: Missing $BAM — skipping $TERA"
            ALL_EXIST=false
            break
        fi
        INPUTS="$INPUTS $BAM"
    done

    if [ "$ALL_EXIST" = false ]; then
        continue
    fi

    OUTBAM="${TERA}_merged.bam"
    if [ -f "$OUTBAM" ]; then
        echo "$OUTBAM already exists, skipping"
        continue
    fi

    echo "=== $(date) === Merging $TERA: ${SRRS[*]} ==="
    "$SAMTOOLS" merge -@ 4 "$OUTBAM" $INPUTS
    "$SAMTOOLS" index -@ 4 "$OUTBAM"
    echo "=== $(date) === Done: $OUTBAM ($(ls -lh "$OUTBAM" | awk '{print $5}')) ==="
done

echo "=== $(date) === All merges complete ==="
