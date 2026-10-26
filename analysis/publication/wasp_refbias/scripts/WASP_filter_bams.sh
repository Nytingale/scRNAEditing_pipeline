#!/bin/bash
# ============================================================================
# WASP_filter_bams.sh — Stage 0b: WASP-tagged BAMs -> WASP-passed BAMs
# ============================================================================
# Publication copy of human/jzlxv3/scripts/wasp_filter_bams.sh, generalised to
# any dataset in WASP_datasets.tsv.
#
# For each WASP-aligned BAM, writes a filtered version keeping only reads that
# either:
#   (a) do not overlap any variant site (no vW tag), OR
#   (b) passed WASP filtering (vW:i:1)
# i.e. it removes the reference-biased reads (vW:i:2-7) while retaining every
# read that never touched an editing site. That retention matters: WASP only
# evaluates reads overlapping a VCF variant, so dropping the untouched reads
# would conflate "bias-corrected" with "restricted to editing sites".
#
# Also writes wasp_filter_summary.csv, the per-unit vW breakdown that Stage 2
# turns into the "loss among WASP-EVALUATED reads" row of WASP_readloss.tsv —
# the number the supplement must quote, not the diluted genome-wide one.
#
# Usage:  qsub -N WASP_filter_<ds> publication/wasp_refbias/scripts/WASP_filter_bams.sh <dataset>
#         bash publication/wasp_refbias/scripts/WASP_filter_bams.sh <dataset>
#
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 4
#$ -l mem=16G
#$ -l h_rt=24:00:00
# Do NOT pass -V (SGE JSV / truncated PATH).

set -euo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

THREADS="${NSLOTS:-4}"

# SGE COPIES THE JOB SCRIPT into a spool dir before running it, so under qsub
# ${BASH_SOURCE[0]} is /cm/local/.../job_scripts/<jobid> and the sibling
# WASP_dataset_config.sh is not there. Prefer the real location when running
# interactively, fall back to the canonical path under qsub.
WASP_SCRIPTS_DEFAULT="${SCRNA_DIR:-.}/publication/wasp_refbias/scripts"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/WASP_dataset_config.sh" ] || SCRIPT_DIR="${WASP_SCRIPTS_DEFAULT}"
[ -f "${SCRIPT_DIR}/WASP_dataset_config.sh" ] || {
    echo "ERROR: cannot locate WASP_dataset_config.sh (tried ${SCRIPT_DIR})" >&2; exit 1; }

# shellcheck source=WASP_dataset_config.sh
source "${SCRIPT_DIR}/WASP_dataset_config.sh"
wasp_load_dataset "${1:-jzlxv3}"

[ -d "$ALLDIR" ] || { echo "ERROR: missing $ALLDIR — run stage 0 first" >&2; exit 1; }
mkdir -p "$PASSDIR"

SUMMARY="${PASSDIR}/wasp_filter_summary.csv"
echo "unit,total_reads,no_variant_overlap,wasp_passed,wasp_failed,pct_removed" > "$SUMMARY"

mapfile -t BAMS < <(find "$ALLDIR" -maxdepth 1 -name '*Aligned.sortedByCoord.out.bam' | sort)
N=${#BAMS[@]}
[ "$N" -ge 1 ] || { echo "ERROR: no WASP BAMs in $ALLDIR" >&2; exit 1; }
echo "=== $(date) === filtering ${N} BAM(s) for ${DATASET} ==="

i=0
for INBAM in "${BAMS[@]}"; do
    i=$((i + 1))
    UNIT=$(basename "$INBAM" _Aligned.sortedByCoord.out.bam)
    OUTBAM="${PASSDIR}/${UNIT}_Aligned.sortedByCoord.out.bam"

    # Idempotent, like stage 0: a 768-unit dataset will be resumed at least
    # once, and re-filtering a finished BAM is pure waste.
    if [ ! -s "$OUTBAM" ] || ! samtools quickcheck "$OUTBAM" 2>/dev/null; then
        samtools view -h -b -@ "$THREADS" -e '![vW] || [vW]==1' "$INBAM" -o "$OUTBAM"
        samtools index -@ "$THREADS" "$OUTBAM"
    fi

    # Totals from the INDEX, not a full pass: idxstats sums mapped+unmapped per
    # reference and was verified to match `samtools view -c` exactly on these
    # BAMs. That removes two whole reads of every file, which matters at 768
    # units even though each is only ~48 MB.
    TOTAL=$(samtools idxstats "$INBAM"  | awk '{s+=$3+$4} END{print s+0}')
    KEPT=$(samtools idxstats "$OUTBAM"  | awk '{s+=$3+$4} END{print s+0}')
    REMOVED=$((TOTAL - KEPT))
    WASP_PASS=$(samtools view -c -@ "$THREADS" -d vW:1 "$INBAM" 2>/dev/null || echo 0)
    NO_VW=$((KEPT - WASP_PASS))

    if [ "$TOTAL" -gt 0 ]; then
        PCT=$(awk "BEGIN {printf \"%.4f\", $REMOVED/$TOTAL * 100}")
    else
        PCT="0"
    fi

    echo "[$i/$N] $UNIT | total $TOTAL | kept $KEPT | removed $REMOVED ($PCT%)"
    echo "${UNIT},${TOTAL},${NO_VW},${WASP_PASS},${REMOVED},${PCT}" >> "$SUMMARY"
done

echo
echo "=== $(date) === Stage 0b DONE (${DATASET}) ==="
echo "Filtered BAMs: $PASSDIR"
echo "Summary CSV:   $SUMMARY"
awk -F, 'NR>1 {t+=$2; p+=$4; f+=$5} END {
    if (t>0) printf "Totals: %d reads, %d evaluated, %d failed (%.2f%% of evaluated)\n",
                    t, p+f, f, (p+f>0 ? 100*f/(p+f) : 0)
}' "$SUMMARY"
echo "Next: qsub publication/wasp_refbias/scripts/WASP_prep_pileups.sh ${DATASET}"
