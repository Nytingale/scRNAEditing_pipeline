#!/bin/bash
# Merge T1–T4 per-arm BAMs into one pseudobulk BAM per CRISPR arm.
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -N merge_gt
#$ -pe threaded 4
#$ -l h_vmem=8G
#$ -l h_rt=8:00:00

set -euo pipefail

export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

BASE="${SCRNA_DIR:-.}/human/pm"
INDIR="${BASE}/02_gt_bams"
ARMS="ADAR1_KO ADARB1_KO ADARB2_KO AAVS1_ctrl NTC_ctrl"

for ARM in $ARMS; do
    OUT="${INDIR}/merged_${ARM}.bam"
    INPUTS=()
    for S in T1 T2 T3 T4; do
        F="${INDIR}/${S}_${ARM}.bam"
        if [ ! -s "$F" ]; then
            echo "ERROR: missing $F — split jobs not done" >&2
            exit 1
        fi
        INPUTS+=("$F")
    done
    if [ -f "$OUT" ] && [ -f "${OUT}.bai" ]; then
        echo "$OUT exists, skipping"
        continue
    fi
    echo "=== $(date) === merging $ARM ==="
    samtools merge -@ 4 -f "$OUT" "${INPUTS[@]}"
    samtools index -@ 4 "$OUT"
    echo "=== $(date) === $OUT  $(samtools view -c "$OUT") reads ==="
done
