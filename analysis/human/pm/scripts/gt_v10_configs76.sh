#!/bin/bash
# gt_v10_configs76.sh   (2026-11-14; was tmp/resp_v10/configs76_seeded.sh)
# A-to-G sites, A-to-G % and mean F1 of the 76 ranked configurations, both datasets (human/scripts/gt_v10_configs76.R)
#   -> human/pm/gt_eval_v9/configs76/configs76_seeded.tsv
# Usage : qsub human/pm/scripts/gt_v10_configs76.sh   (from the repo root)
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -o ${SCRNA_DIR:-.}/human/pm/gt_eval_v9/configs76/
#$ -j y
cd ${SCRNA_DIR:-.}
${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript human/scripts/gt_v10_configs76.R
