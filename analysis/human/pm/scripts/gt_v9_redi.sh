#!/bin/bash
# gt_v9_redi.sh STAGE=universe|ucov|score|splice|submit   (user 2026-11-07)
# REDIportal truths under the v9 rules (human/scripts/gt_score_rediportal_v9.R), state mismatch, 60 configurations,
# pwcov settings -> human/pm/gt_eval_v9/rediportal_v9/
#   universe  GT_REDI_PART=universe -> R/universe/; unit_cov setup at the reverse-orientation positives
#   ucov      array 1-28: gt_v8_unit_cov.sh LIB=reverse -> R/unit_cov_rev/ko_truth/<bam>.tsv.gz
#   score     array 1-2: published | Cellsnp-lite hybrid -> R/_parts/{pub,ch}/scoring_mismatch/
#   splice    gt_v9_redi_splice.R -> R/scoring_mismatch/
#   submit    login node (qsub only)
#$ -S /bin/bash
#$ -N gt_v9_redi
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/human/pm/gt_eval_v9/logs/
#$ -e ${SCRNA_DIR:-.}/human/pm/gt_eval_v9/logs/
#$ -l h_vmem=48G
#$ -l h_rt=12:00:00
set -euo pipefail
export LC_ALL=C
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
export GT_V9_DIR=$PM/gt_eval_v9
R=$GT_V9_DIR/rediportal_v9
SELF=$PM/scripts/gt_v9_redi.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
export GT_REDI_UDIR=$R/universe GT_REDI_REVCOV=$R/unit_cov_rev GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pwcov.txt
export GT_REDI_EXCLUDE_BED=$PM/pgp1_snps/PGP1_hu43860C_hg38_snvs.bed GT_CALL_RE='_step5_nomm\.tsv\.gz$'
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
STAGE=${STAGE:-${1:-}}
cd "$B"

case $STAGE in
universe)
  GT_REDI_PART=universe $RS human/scripts/gt_score_rediportal_v9.R
  OUT=$R/unit_cov_rev TRUTH=$R/universe/rev_positive_sites.tsv.gz bash $PM/scripts/gt_v8_unit_cov.sh setup
  ;;
ucov) OUT=$R/unit_cov_rev LIB=reverse bash $PM/scripts/gt_v8_unit_cov.sh ;;
score)
  if [ "$SGE_TASK_ID" = 1 ]; then d="arm=${PUB[arm]},sample=${PUB[sample]},celltype=${PUB[celltype]}"; o=$R/_parts/pub/scoring_mismatch
  else c=gt_eval_v8/cellsnp_hybrid/results/centre; o=$R/_parts/ch/scoring_mismatch
    d="arm=$c/${PUB[arm]},sample=$c/${PUB[sample]},celltype=$c/${PUB[celltype]}"; fi
  GT_REDI_PART=score GT_RES_DIRS="$d" GT_REDI_OUT=$o $RS human/scripts/gt_score_rediportal_v9.R
  ;;
splice) $RS human/scripts/gt_v9_redi_splice.R ;;
submit)
  mkdir -p $GT_V9_DIR/logs
  Q=$(command -v qsub || echo /cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub)
  j1=$($Q -terse -v STAGE=universe $SELF)
  j2=$($Q -terse -v STAGE=ucov -t 1-28 -l h_vmem=4G -hold_jid $j1 $SELF | cut -d. -f1)
  j3=$($Q -terse -v STAGE=score -t 1-2 -l h_vmem=32G -hold_jid $j2 $SELF | cut -d. -f1)
  j4=$($Q -terse -v STAGE=splice -l h_vmem=8G -hold_jid $j3 $SELF)
  echo "universe $j1 ucov $j2 score $j3 splice $j4"
  ;;
*) echo "STAGE?" >&2; exit 1 ;;
esac
