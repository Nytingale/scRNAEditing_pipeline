#!/bin/bash
# mding_3prime_bias_count.sh   (2026-11-07) mding_3prime_bias.R per mouse x Smart-seq2 seed
#   tasks 1-6: mouse (Cortex1: 1-3, Cortex2: 4-6) x seed
#   qsub -t 1-6 -hold_jid mding_pus mouse/scripts/mding_3prime_bias_count.sh
#$ -S /bin/bash
#$ -N mding_cnt
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/mjzl10x/three_prime_bias/logs/
#$ -j y
#$ -l h_vmem=32G
#$ -l h_rt=6:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
MICE=(Cortex1 Cortex2); SEEDS=(20261024 20261025 20261026)
i=$((SGE_TASK_ID - 1))
MOUSE=${MICE[$((i / 3))]} SEED=${SEEDS[$((i % 3))]} ${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript $B/mouse/scripts/mding_3prime_bias.R
