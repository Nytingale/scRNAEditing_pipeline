#!/bin/bash
# gt_v9_build.sh STAGE=truth|ucov|collect|score|rank|rank76|submit   (user 2026-11-07)
# Teratoma ground truth v9 (human/scripts/gt_truth_v9.R) and the headline rescore (state mismatch, 60
# configurations, three tiers, pwcov settings, GT_FP_RULE=ag) into human/pm/gt_eval_v9/.
#   truth    gt_truth_v7.R + gt_truth_v7_variants.R (GT_V7_VARIANT=AGTC, w1 searched) -> V9/build/;
#            gt_truth_v9.R -> V9/truth/; unit_cov setup (positives only)
#   ucov     array 1-56: gt_v8_unit_cov.sh, tasks 1-28 forward (read A>G), 29-56 reverse (read T>C)
#   collect  gt_v9_unit_cov_collect.R -> V9/unit_cov/ko_truth_unit_cov.tsv.gz
#   score    array 1-6: tier x (58 published | 2 Cellsnp-lite hybrid) -> V9/_parts/{pub,ch}/mismatch/<tier>/
#   rank     gt_v9_splice_compare.R splice, gt_rank_v7.R, gt_v9_splice_compare.R compare
#   rank76   after umimd_score.sh: gt_v9_splice_compare.R splice76, gt_rank_v7.R -> V9/scoring_pwcov76/mismatch/
#   submit   login node (qsub only)
#$ -S /bin/bash
#$ -N gt_v9
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
C=$V8/cellsnp_hybrid
export GT_V9_DIR=${GT_V9_DIR:-$PM/gt_eval_v9}
V9=$GT_V9_DIR
SELF=$PM/scripts/gt_v9_build.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
TIERS=(arm sample celltype)
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
export GT_UNIT_RULE=gate_cov GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pwcov.txt GT_FP_RULE=ag
STAGE=${STAGE:-${1:-}}
cd "$B"

case $STAGE in
truth)
  GT_V7_VARIANT=AGTC GT_V7_TRUTH_DIR=$V9/build/v7truth $RS human/scripts/gt_truth_v7.R
  GT_V7_VARIANT=AGTC GT_V7_TRUTH_DIR=$V9/build/v7truth GT_V7_VARIANTS_DIR=$V9/build/variants \
    $RS human/scripts/gt_truth_v7_variants.R
  $RS human/scripts/gt_truth_v9.R
  for l in fwd rev; do OUT=$V9/unit_cov_$l TRUTH=$V9/truth/gt_v9_positive_sites.tsv.gz bash $PM/scripts/gt_v8_unit_cov.sh setup; done
  ;;
ucov)
  k=$(( SGE_TASK_ID - 1 ))
  if [ $k -lt 28 ]; then l=fwd; lib=forward; else l=rev; lib=reverse; fi
  SGE_TASK_ID=$(( k % 28 + 1 )) OUT=$V9/unit_cov_$l LIB=$lib bash $PM/scripts/gt_v8_unit_cov.sh
  ;;
collect)
  for l in fwd rev; do OUT=$V9/unit_cov_$l bash $PM/scripts/gt_v8_unit_cov.sh collect; done
  $RS human/scripts/gt_v9_unit_cov_collect.R
  ;;
score)
  i=$(( SGE_TASK_ID - 1 )); tier=${TIERS[$(( i % 3 ))]}
  if [ $i -lt 3 ]; then res=$PM/${PUB[$tier]}; o=$V9/_parts/pub/mismatch/$tier
  else res=$C/results/centre/${PUB[$tier]}; o=$V9/_parts/ch/mismatch/$tier; fi
  GT_V7_TRUTH=$V9/truth/gt_v9_truth_sites.tsv.gz GT_UNIT_COV=$V9/unit_cov/ko_truth_unit_cov.tsv.gz \
    $RS human/scripts/gt_score_v7.R "$tier" "$res" "$o" '_step5_nomm\.tsv\.gz$'
  ;;
rank)
  $RS human/scripts/gt_v9_splice_compare.R splice
  $RS human/scripts/gt_rank_v7.R $V9/scoring_pwcov/mismatch > $V9/scoring_pwcov/mismatch/rank.log 2>&1
  $RS human/scripts/gt_v9_splice_compare.R compare
  ;;
rank76)
  $RS human/scripts/gt_v9_splice_compare.R splice76
  $RS human/scripts/gt_rank_v7.R $V9/scoring_pwcov76/mismatch > $V9/scoring_pwcov76/mismatch/rank.log 2>&1
  ;;
submit)
  mkdir -p $V9/logs
  Q=$(command -v qsub || echo /cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub)
  j1=$($Q -terse -v STAGE=truth -l h_vmem=32G $SELF)
  j2=$($Q -terse -v STAGE=ucov -t 1-56 -l h_vmem=4G -hold_jid $j1 $SELF | cut -d. -f1)
  j3=$($Q -terse -v STAGE=collect -l h_vmem=16G -hold_jid $j2 $SELF)
  j4=$($Q -terse -v STAGE=score -t 1-6 -hold_jid $j3 $SELF | cut -d. -f1)
  j5=$($Q -terse -v STAGE=rank -l h_vmem=16G -hold_jid $j4 $SELF)
  echo "truth $j1 ucov $j2 collect $j3 score $j4 rank $j5"
  ;;
*) echo "STAGE?" >&2; exit 1 ;;
esac
