#!/bin/bash
# Submit annotation only. CellRanger + dual ref must already be done.
# Optional: qsub -hold_jid 1619564 this script's job if CellRanger is still running.

set -euo pipefail

PMDIR="${SCRNA_DIR:-.}/human/pm"
SCRIPTS="$PMDIR/scripts"
LOGS="$PMDIR/logs"
mkdir -p "$LOGS"

if [ "${1:-}" != "" ]; then
  echo "Will hold on job $1"
  JOB=$(qsub -hold_jid "$1" "$SCRIPTS/run_teratoma_annotation.sh" | awk '{print $3}')
else
  JOB=$(qsub "$SCRIPTS/run_teratoma_annotation.sh" | awk '{print $3}')
fi

echo "Submitted ter_annot: $JOB"
echo "Log: $LOGS/teratoma_annotation.log"
