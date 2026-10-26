#!/bin/bash
# depmap_adar_isg_exADAR.sh   (2026-10-04; Reviewer 3, Major 10)
# As depmap_adar_isg.sh (DepMap 26Q1 lung lines), but the ISG scores exclude only ADAR (the knocked-out
# gene); dsRNA sensors stay in the gene sets. Then the lungTME figure in which panel c ranks single ISGs
# and sensors by their correlation with the ADAR KO effect.
# Output: human/dl/depmap_results_exADAR/, publication/lungTME/figures/dl/LUNGTME_panel_depmap_adar_isg_exADAR.svg
#$ -S /bin/bash
#$ -N depmap_exADAR
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=24G
#$ -l h_rt=4:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
R=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
G=$B/human/dl/genesets
[ -s $G/ISGs.tsv ] || { echo "missing $G/ISGs.tsv" >&2; exit 1; }
export ISG_FILES="$G/ds_hallmark_sets.gmt;HCNetlas_423_ISGs=$G/ISGs.tsv"
export ISG_EXCLUDE=ADAR
export DEPMAP_OUT=$B/human/dl/depmap_results_exADAR
$R $B/human/scripts/depmap_adar_isg.R
$R $B/publication/lungTME/scripts/LUNGTME_depmap_plot_exADAR.R
