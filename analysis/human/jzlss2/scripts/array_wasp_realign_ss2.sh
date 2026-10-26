#!/bin/bash
# ============================================================================
# Array task wrapper for Smart-seq2 STAR-WASP re-alignment
# Called by SGE array job — $SGE_TASK_ID indexes into srr_array_list.txt
# ============================================================================
# Same list and same alignment script as array_align_ss2.sh; the only
# difference is WASP_VCF + OUTDIR, which switch align_ss2.sh into WASP mode.
# Keeping one alignment script means the "all reads" and "WASP-passed" arms of
# the comparison cannot differ in anything except the WASP flags.
# ============================================================================

set -euo pipefail

# Activate the env here rather than relying on qsub -V carrying it in: a 768-
# task array that inherits a half-set-up shell fails 768 times. Same preamble
# as WASP_prep_pileups.sh — PS1 and the `|| true` are needed because the
# interactive rc touches unbound variables and would abort under `set -u`.
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:${PATH:-}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

command -v STAR     >/dev/null || { echo "ERROR: STAR not on PATH" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }

BASEDIR="${SCRNA_DIR:-.}/human/jzlss2"
REFDIR="${SCRNA_DIR:-.}/human/ref_annots"

export WASP_VCF="${REFDIR}/REDIportal_hg38_v2.vcf"
export OUTDIR="${BASEDIR}/01_starsolobams_wasp"
export NTHREADS="${NSLOTS:-2}"

SRR=$(sed -n "${SGE_TASK_ID}p" "${BASEDIR}/srr_array_list.txt")

if [ -z "$SRR" ]; then
    echo "ERROR: No SRR found for task ID $SGE_TASK_ID"
    exit 1
fi

echo "Task $SGE_TASK_ID: WASP re-aligning $SRR"
bash "${BASEDIR}/scripts/align_ss2.sh" "$SRR"
