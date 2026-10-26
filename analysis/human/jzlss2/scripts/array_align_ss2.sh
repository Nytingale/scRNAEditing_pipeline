#!/bin/bash
# ============================================================================
# Array task wrapper for Smart-seq2 STAR alignment
# Called by SGE array job — $SGE_TASK_ID indexes into srr_array_list.txt
# ============================================================================

set -euo pipefail

BASEDIR="${SCRNA_DIR:-.}/human/jzlss2"
SRR=$(sed -n "${SGE_TASK_ID}p" "${BASEDIR}/srr_array_list.txt")

if [ -z "$SRR" ]; then
    echo "ERROR: No SRR found for task ID $SGE_TASK_ID"
    exit 1
fi

echo "Task $SGE_TASK_ID: aligning $SRR"
bash "${BASEDIR}/scripts/align_ss2.sh" "$SRR"
