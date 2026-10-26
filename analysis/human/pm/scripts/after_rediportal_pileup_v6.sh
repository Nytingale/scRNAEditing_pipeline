#!/bin/bash
# after_rediportal_pileup_v6.sh   (HANDOFF 2026-09-29 UPDATE-25)
# After rediportal_pileup_v6.sh: check the merged-NTC pileup, flag REDIportal positions in
# dbSNP156 common (same BED and bedtools as the v6 candidates), score against REDIportal,
# draw the comparison figures and the cross-tier grid under each REDIportal truth variant.
#   qsub -hold_jid redi_pileup after_rediportal_pileup_v6.sh
#$ -S /bin/bash
#$ -N after_redi
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
set -euo pipefail
B=${SCRNA_DIR:-.}
V6=$B/human/pm/gt_eval_v6
RD=$V6/rediportal
R=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
S=$B/publication/teratoma_gt/scripts
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
n=$(find "$RD/pileup" -name 'chr*.tsv.gz' | wc -l)
[ "$n" -eq "$(wc -l < "$RD/tasks.tsv")" ] || { echo "ERROR: $n pileup outputs, expected $(wc -l < "$RD/tasks.tsv")" >&2; exit 1; }
echo "bedtools: $(bedtools --version)"
cat "$RD"/bed/chr*.bed | sort -k1,1 -k2,2n -S 8G \
  | bedtools intersect -u -a - -b "$V6/snp/dbSNP156_common.bed" \
  | awk -v OFS='\t' '{print $1 ":" $3}' > "$RD/rediportal_dbsnp_common_sites.txt"
echo "REDIportal positions in dbSNP156 common: $(wc -l < "$RD/rediportal_dbsnp_common_sites.txt")"
$R $B/human/scripts/gt_score_rediportal_v6.R
$R $S/TERATOMA_rediportal_compare_plot.R
declare -A TITLE=([cov10_alt1]="REDIportal sites, >= 10 reads, >= 1 G in pooled controls" \
                  [alt3]="REDIportal sites, >= 10 reads, >= 3 G in pooled controls" \
                  [alt3_e05]="REDIportal sites, >= 10 reads, >= 3 G, edit >= 5% in pooled controls")
for v in cov10_alt1 alt3 alt3_e05; do
  GRID_SCORING_TSV=$RD/scoring/gt_tier_scoring_$v.tsv GRID_STEM=TERATOMA_cross_tier_grid_rediportal_$v \
  GRID_TITLE_SUFFIX=" (truth: ${TITLE[$v]})" $R $S/TERATOMA_cross_tier_grid_plot.R
done
$R $S/TERATOMA_rank_robustness_plot.R
