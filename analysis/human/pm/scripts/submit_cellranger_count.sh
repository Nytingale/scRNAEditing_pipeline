#!/bin/bash
# Submit CellRanger count array for corrected T1–T4 pairings.
# Reference already built; this only remakes count.

set -euo pipefail

PMDIR="${SCRNA_DIR:-.}/human/pm"
SCRIPTS="$PMDIR/scripts"
LOGS="$PMDIR/logs"
mkdir -p "$LOGS"

JOB_CR=$(qsub -N cr_pm -t 1-4 -pe threaded 16 -l h_vmem=8G -l h_rt=24:00:00 \
  -j y -o "$LOGS/cellranger_count_\$TASK_ID.log" -cwd -V -b y \
  "bash $SCRIPTS/run_cellranger_count.sh" | awk '{print $3}')
echo "Submitted cr_pm (array 1-4): $JOB_CR"
echo "Monitor with: qstat -u \$USER"
