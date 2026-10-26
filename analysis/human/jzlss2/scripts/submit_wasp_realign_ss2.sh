#!/bin/bash
# ============================================================================
# Submit the STAR-WASP re-alignment array for Smart-seq2 (JZLSS2)
# ============================================================================
# Stage 0 of the WASP reference-bias diagnostic on a SECOND dataset
# (Reviewer 2, Major 4). 768 SRRs, one cell each.
#
# Usage:
#   bash submit_wasp_realign_ss2.sh pilot        # tasks 1-8 only, verify first
#   bash submit_wasp_realign_ss2.sh              # the full 1-768
#   TC=12 bash submit_wasp_realign_ss2.sh        # override the throttle
#
# >>> RUN THE PILOT FIRST. <<<
# The plain SS2 alignment already peaked at 30.89 GB of a 32.00 GB cap (96.5%,
# job 1599415.100), and STAR-WASP additionally holds the REDIportal VCF
# (~16 M sites) in memory. So the existing `-pe threaded 2 -l h_vmem=16G`
# WOULD OOM on every one of the 768 tasks. The request below is sized up;
# confirm MEM_USED_G against MEM_LIMIT_G in a pilot log before spending 768
# index loads on an assumption.
#
# THROTTLE. -tc caps concurrent tasks. The binding constraint is not slots but
# I/O: every task loads the 31 GB hg38 index from shared scratch, so 768
# unthrottled tasks are a filesystem storm, not just a queue full. The GT
# Nextflow runs also need slots, and they submit continuously. TC=20 at 2
# slots each = 40 slots, which leaves the cluster usable.
# ============================================================================

set -euo pipefail

BASEDIR="${SCRNA_DIR:-.}/human/jzlss2"
LOGDIR="${BASEDIR}/logs"
LIST="${BASEDIR}/srr_array_list.txt"
mkdir -p "$LOGDIR"

MODE="${1:-full}"
TC="${TC:-20}"

[ -s "$LIST" ] || { echo "ERROR: missing $LIST" >&2; exit 1; }
N_SRRS=$(wc -l < "$LIST")

case "$MODE" in
    pilot) RANGE="1-8";        NAME="wasp_ss2_pilot" ;;
    full)  RANGE="1-${N_SRRS}"; NAME="wasp_ss2" ;;
    *)     RANGE="$MODE";       NAME="wasp_ss2" ;;   # e.g. "9-768"
esac

echo "SRRs in list : ${N_SRRS}"
echo "task range   : ${RANGE}"
echo "throttle     : ${TC} concurrent"

# Memory: `-l mem` is the knob the cluster actually enforces (see CHANGELOG,
# "SCHEDULER MEMORY REQUESTS"); h_vmem is kept because it is per-slot here and
# 2 x 24G is what produced the observed 32 G cap on the plain run. Both are set
# to the same 48 G so declared == enforced either way.
# h_rt: the plain alignment took 5.4 min/cell; WASP re-maps every
# variant-overlapping read, so 4 h is generous headroom, not an estimate.
qsub -N "$NAME" \
    -t "$RANGE" \
    -tc "$TC" \
    -pe threaded 2 \
    -l h_rt=4:00:00 \
    -l h_vmem=24G \
    -l mem=48G \
    -j y \
    -o "${LOGDIR}/${NAME}_\$TASK_ID.log" \
    -cwd -V \
    "${BASEDIR}/scripts/array_wasp_realign_ss2.sh"

cat <<EOF

Submitted ${NAME}, tasks ${RANGE}, throttled to ${TC} concurrent.

Verify BEFORE releasing the rest (if this was the pilot):
  grep -E 'MEM_LIMIT_G|MEM_USED_G|ELAPSED_TIME_S|JOB_EXIT_STATUS' ${LOGDIR}/wasp_ss2_pilot_1.log
  # MEM_USED_G must sit clearly below MEM_LIMIT_G, and exit status 0.

Then:
  bash ${BASEDIR}/scripts/submit_wasp_realign_ss2.sh 9-${N_SRRS}

Progress:
  ls ${BASEDIR}/01_starsolobams_wasp/*Aligned.sortedByCoord.out.bam | wc -l   # target ${N_SRRS}

Re-running is safe: align_ss2.sh skips any cell whose BAM already passes
samtools quickcheck, so a resubmission only picks up the gaps.
EOF
