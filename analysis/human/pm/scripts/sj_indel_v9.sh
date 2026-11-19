#!/bin/bash
# Response Figure 17 on the v9 truth: KO ground truth sites (t = 50%) by distance to the nearest splice junction and
# indel. Indels re-extracted from the seed-8787 merged NTC BAM (ntc_indel_positions_v8.sh), then
# gt_v8_truth_sj_indel_nt.R and GTDEF_v8_sj_indel_{source_data,plot}.R with the v9 paths.
#   STAGE=submit   (login node: qsub only)     STAGE=summary   (held on the indel array)
#$ -N sjind9
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
set -eo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
OUT=$PM/gt_eval_v9/indel_sj
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
case "${STAGE:?}" in
submit)
  mkdir -p "$OUT/indels"
  j1=$(qsub -terse -t 1-23 -v BAM="$PM/02_gt_bams_umi_v6/merged_NTC_ctrl.bam",INDEL_OUT="$OUT/indels" \
       human/pm/scripts/ntc_indel_positions_v8.sh | cut -d. -f1)
  j2=$(qsub -terse -hold_jid "$j1" -v STAGE=summary human/pm/scripts/sj_indel_v9.sh)
  echo "indels $j1 summary $j2"
  ;;
summary)
  GT_TRUTH=$PM/gt_eval_v9/truth/gt_v9_truth_sites.tsv.gz SJ_INDEL_OUT=$OUT $RS human/scripts/gt_v8_truth_sj_indel_nt.R
  S=publication/gt_definition/scripts
  SJ_INDEL_DIR=$OUT $RS $S/GTDEF_v8_sj_indel_source_data.R
  $RS $S/GTDEF_v8_sj_indel_plot.R
  ;;
esac
