#!/bin/bash
# =============================================================================
# gt_robustness_v9.sh   (2026-11-16)
# v9 ground-truth robustness (Response Figure 1) (centre neighbours, m and KO-coverage curves): human/scripts/gt_robustness_v9.R
# Outputs: human/pm/gt_eval_v9/robustness/   Logs: logs/
# Usage  : qsub human/pm/scripts/gt_robustness_v9.sh
#$ -cwd
#$ -N gt_robust_v9

#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=48G
#$ -l h_rt=4:00:00

set -eo pipefail
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true
B="${SCRNA_DIR:-.}"
cd "$B"
/usr/bin/time -v ${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript human/scripts/gt_robustness_v9.R
