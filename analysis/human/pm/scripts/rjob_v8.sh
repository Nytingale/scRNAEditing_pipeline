#!/bin/bash
# rjob_v8.sh <script.R | script.sh> [args...]   (2026-10-15)
# Runs one R or bash analysis step as an SGE job from the project root (Rscript of the scRNAEditing env).
#   qsub [-pe threaded N] [-l h_vmem=..] [-hold_jid ..] -N name human/pm/scripts/rjob_v8.sh human/scripts/x.R args
#$ -S /bin/bash
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
set -euo pipefail
cd ${SCRNA_DIR:-.}
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
echo "=== $(date) $*"
case "$1" in
  *.R) Rscript "$@" ;;
  *)   bash "$@" ;;
esac
echo "=== $(date) done"
