#!/bin/bash
# =============================================================================
# aorta_score_v8.sh  STAGE=submit|score|rank|compare   (2026-10-26)
#
# Score the aorta caller runs (run_pipeline_aorta_v8.sh -> results_v8/) against the aorta truth
# (truth/eval/excluded/gt_aorta_truth_sites.tsv.gz, MGP SNP positions removed; REDIportal negatives
# unscored) as the teratoma headline (gt_eval_v8/scoring_pwcov): gt_score_v7.R with GT_UNIT_RULE=gate_cov
# (per-unit recall denominator: power-gated positives need unit coverage x merged CTRL editing >= 3 and
# >= 10 reads; depth-gated >= 10 reads), unit coverage from truth/eval/excluded/unit_cov.tsv.gz.
#   tiers   arm = CTRL (merged 3 control mice), sample = the 3 control mice (no cell-type tier)
#   states  mismatch = *_step5_nomm (headline), none = *_step4_cov_edit_filtered
#   rank    gt_rank_v7.R with GT_RANK_TIERS=arm,sample;  compare: aorta_rank_compare.R
# Outputs: gt_eval_v8/scoring_pwcov/<state>/{arm,sample}/, config_rank.tsv, caller_rank.tsv; gt_eval_v8/compare/
#   STAGE=submit bash scripts/aorta_score_v8.sh [hold_jid]   (login node: collect -> score 1-4 -> rank -> compare)
#   Env: TRUTH_SET (truth/eval/<set>, default excluded), SC_NAME (gt_eval_v8/<name>, default scoring_pwcov),
#        CMP_NAME (default compare), REPLOGLE_TRUTH (default v7_nomm_k562), NOCOLLECT=1 skips collect,
#        GT_DROP_SITES (sites whose calls gt_score_v7.R removes)
# =============================================================================
#$ -N aorta_score
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
set -eo pipefail
B="${SCRNA_DIR:-.}"
A="$B/mouse/aorta_gse255661"
SELF="$A/scripts/aorta_score_v8.sh"
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
RES="$A/results_v8"
TRUTH_SET=${TRUTH_SET:-excluded}; SC_NAME=${SC_NAME:-scoring_pwcov}; CMP_NAME=${CMP_NAME:-compare}
REPLOGLE_TRUTH=${REPLOGLE_TRUTH:-v7_nomm_k562}
TR="$A/truth/eval/$TRUTH_SET"
SC="$A/gt_eval_v8/$SC_NAME"
STATES=(mismatch none); TIERS=(arm sample)
STAGE=${STAGE:?STAGE=submit|score|rank|compare}
cd "$B"

if [ "$STAGE" = "submit" ]; then
  H=${1:+-hold_jid $1}
  V="TRUTH_SET=$TRUTH_SET,SC_NAME=$SC_NAME,CMP_NAME=$CMP_NAME,REPLOGLE_TRUTH=$REPLOGLE_TRUTH,GT_DROP_SITES=${GT_DROP_SITES:-}"
  if [ "${NOCOLLECT:-0}" = 1 ]; then jc=none; HS=$H
  else jc=$(qsub -terse -N aorta_collect -l h_vmem=8G $H -v STAGE=collect "$A/scripts/run_pipeline_aorta_v8.sh"); HS="-hold_jid $jc"; fi
  js=$(qsub -terse -t 1-4 $HS -v STAGE=score,$V "$SELF" | cut -d. -f1)
  jr=$(qsub -terse -N aorta_rank -l h_vmem=8G -hold_jid "$js" -v STAGE=rank,$V "$SELF")
  jm=$(qsub -terse -N aorta_cmp -l h_vmem=8G -hold_jid "$jr" -v STAGE=compare,$V "$SELF")
  echo "collect $jc | score $js (1-4) | rank $jr | compare $jm"
  exit 0
fi

if [ "$STAGE" = "score" ]; then
  i=$(( SGE_TASK_ID - 1 )); st=${STATES[$(( i / 2 ))]}; tier=${TIERS[$(( i % 2 ))]}
  re=$( [ "$st" = mismatch ] && echo '_step5_nomm\.tsv\.gz$' || echo '_step4_cov_edit_filtered\.tsv\.gz$' )
  ure=$( [ "$tier" = arm ] && echo '^CTRL$' || echo '^Control_' )
  echo "state=$st tier=$tier units=$ure"
  GT_V7_TRUTH="$TR/gt_aorta_truth_sites.tsv.gz" GT_UNIT_COV="$TR/unit_cov.tsv.gz" GT_UNIT_RULE=gate_cov \
  GT_UNIT_EDIT_UNIT=CTRL GT_UNIT_RE="$ure" \
    $RS human/scripts/gt_score_v7.R "$tier" "$RES" "$SC/$st/$tier" "$re"
  exit 0
fi

if [ "$STAGE" = "rank" ]; then
  for st in "${STATES[@]}"; do GT_RANK_TIERS=arm,sample $RS human/scripts/gt_rank_v7.R "$SC/$st"; done
  exit 0
fi

if [ "$STAGE" = "compare" ]; then
  AORTA_SCORE_DIR="$SC_NAME" AORTA_CMP_DIR="$CMP_NAME" REPLOGLE_TRUTH="$REPLOGLE_TRUTH" $RS "$A/scripts/aorta_rank_compare.R"
  exit 0
fi
echo "ERROR: unknown STAGE '$STAGE'" >&2; exit 1
