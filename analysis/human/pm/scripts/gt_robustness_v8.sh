#!/bin/bash
# =============================================================================
# gt_robustness_v8.sh   (2026-10-13)
# v8 ground-truth robustness (centre neighbours, m and KO-coverage curves): human/scripts/gt_robustness_v8.R
# Outputs: human/pm/gt_eval_v8/robustness/   Logs: logs/
# Usage  : qsub human/pm/scripts/gt_robustness_v8.sh
#$ -cwd
#$ -N gt_robust_v8
#$ -q slice
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=16G
#$ -l h_rt=4:00:00

set -eo pipefail
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true
B="${SCRNA_DIR:-.}"
cd "$B"
/usr/bin/time -v ${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript human/scripts/gt_robustness_v8.R
