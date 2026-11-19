#!/bin/bash
# aorta_v9_build.sh STAGE=truth|score|rank|submit   (user 2026-11-07)
# Aorta ground truth v9 (aorta_truth_v9.R) and the headline rescore (state mismatch, 58 published + 2 Cellsnp-lite
# hybrid configurations, arm and sample tiers, aorta_score_v8.sh settings, GT_FP_RULE=ag) into gt_eval_v9/.
#   truth   aorta_truth_v9.R -> truth/eval_v9/excluded/
#   score   array 1-4: tier x (published | hybrid) -> gt_eval_v9/_parts/{pub,ch}/mismatch/<tier>/
#   rank    aorta_v9_splice_compare.R splice, gt_rank_v7.R (v9 and v8 at 60), compare
#   rank76  after aorta_umimd_score.sh: splice76, gt_rank_v7.R -> V9/scoring_pwcov76/mismatch/ (2026-11-14)
#   submit  qsub only
#$ -S /bin/bash
#$ -N aorta_v9
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/gt_eval_v9/logs/
#$ -j y
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
A=$B/mouse/aorta_gse255661
V9=${AORTA_GT_V9_DIR:-$A/gt_eval_v9}
TR=${AORTA_V9_DIR:-$A/truth/eval_v9/excluded}
SELF=$A/scripts/aorta_v9_build.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
export GT_FP_RULE=ag
TIERS=(arm sample)
STAGE=${STAGE:-${1:-}}
cd "$B"

case $STAGE in
truth) $RS $A/scripts/aorta_truth_v9.R ;;
score)
  i=$(( SGE_TASK_ID - 1 )); tier=${TIERS[$(( i % 2 ))]}
  ure=$( [ "$tier" = arm ] && echo '^CTRL$' || echo '^Control_' )
  if [ $i -lt 2 ]; then res=$A/results_v8; o=$V9/_parts/pub/mismatch/$tier
  else res=$A/gt_eval_v8/cellsnp_hybrid/results/results_v8; o=$V9/_parts/ch/mismatch/$tier; fi
  GT_V7_TRUTH=$TR/gt_aorta_truth_sites.tsv.gz GT_UNIT_COV=$TR/unit_cov.tsv.gz GT_UNIT_RULE=gate_cov \
  GT_UNIT_EDIT_UNIT=CTRL GT_UNIT_RE="$ure" \
    $RS human/scripts/gt_score_v7.R "$tier" "$res" "$o" '_step5_nomm\.tsv\.gz$'
  ;;
rank)
  $RS $A/scripts/aorta_v9_splice_compare.R splice
  for d in $V9/scoring_pwcov/mismatch $V9/v8_60/scoring_pwcov/mismatch; do
    GT_RANK_TIERS=arm,sample $RS human/scripts/gt_rank_v7.R $d > $d/rank.log 2>&1
  done
  $RS $A/scripts/aorta_v9_splice_compare.R compare
  ;;
rank76)
  $RS $A/scripts/aorta_v9_splice_compare.R splice76
  GT_RANK_TIERS=arm,sample $RS human/scripts/gt_rank_v7.R $V9/scoring_pwcov76/mismatch > $V9/scoring_pwcov76/mismatch/rank.log 2>&1
  ;;
submit)
  mkdir -p $V9/logs
  Q=/cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub
  j1=$($Q -terse -v STAGE=score -t 1-4 $SELF | cut -d. -f1)
  j2=$($Q -terse -v STAGE=rank -l h_vmem=16G -hold_jid $j1 $SELF)
  echo "score $j1 rank $j2"
  ;;
*) echo "STAGE?" >&2; exit 1 ;;
esac
