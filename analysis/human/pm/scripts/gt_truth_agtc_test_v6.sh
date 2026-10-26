#!/bin/bash
# =============================================================================
# gt_truth_agtc_test_v6.sh   (2026-10-06)
# Strand-agnostic ground-truth test (A>G vs A>G + T>C vs pooled orientations, multiallelic
# kept / removed) re-scored on the published call sets: human/scripts/gt_truth_agtc_test_v6.R
# Outputs: human/pm/gt_eval_v6/agtc_test/   Logs: logs/
# Usage  : qsub human/pm/scripts/gt_truth_agtc_test_v6.sh
#$ -cwd
#$ -N gt_agtc_test
#$ -q slice
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=12G
#$ -l h_rt=12:00:00

set -eo pipefail
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true
B="${SCRNA_DIR:-.}"
cd "$B"
${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript human/scripts/gt_truth_agtc_test_v6.R
