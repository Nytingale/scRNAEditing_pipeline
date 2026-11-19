#!/bin/bash
# ding_3prime_bias_sub_count.sh   (2026-10-24; Reviewer 3, major 2)
# REDIportal site counting and gene-body positions for the subsampled 10x UMI-aware MarkDuplicates pileup vs each
# subsampled Smart-seq2 coordinate-deduplicated pileup (ding_3prime_bias_pileup_sub.sh), 217,323,582 reads each, no
# thinning (until 2026-11-14: X10=x10us20261024 SS2=ss2p1ds<seed>, umi_tools, 217,323,670 reads).
#   qsub -t 1-3 -hold_jid ding_pus human/scripts/ding_3prime_bias_sub_count.sh
#$ -S /bin/bash
#$ -N ding_subc
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=64G
#$ -l h_rt=6:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
SEEDS=(20261024 20261025 20261026)
export X10=x10ms20261024 SS2=ss2p1dsm${SEEDS[$((SGE_TASK_ID - 1))]} DEPTH=subsampled
R=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
$R $B/human/scripts/ding_3prime_bias_rediportal_pu.R
$R $B/human/scripts/ding_3prime_bias_genebody.R
