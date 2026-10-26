#!/bin/bash
# depmap_adar_isg.sh   (2026-10-04; Reviewer 3, Major 10)
# Runs human/scripts/depmap_adar_isg.R (DepMap 26Q1 ADAR dependency vs interferon state), then the
# lungTME DepMap figure. Output: human/dl/depmap_results/, publication/lungTME/figures/dl/
#$ -S /bin/bash
#$ -N depmap_adar
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=24G
#$ -l h_rt=4:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
R=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
# gene sets of the Down syndrome analysis: Hallmark sets extracted from 13_ISGs.R, HCNetlas 423 ISGs
G=$B/human/dl/genesets
[ -s $G/ISGs.tsv ] || { echo "missing $G/ISGs.tsv" >&2; exit 1; }
export ISG_FILES="$G/ds_hallmark_sets.gmt;HCNetlas_423_ISGs=$G/ISGs.tsv"
$R $B/human/scripts/depmap_adar_isg.R
$R $B/publication/lungTME/scripts/LUNGTME_depmap_plot.R
