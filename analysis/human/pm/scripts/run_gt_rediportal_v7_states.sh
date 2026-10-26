#!/bin/bash
# run_gt_rediportal_v7_states.sh <none|mismatch|clustered>   (2026-10-11)
# REDIportal scoring of the three non-published call-set window states, under gt_eval_v7
# (the published state is run_gt_rediportal_v7.sh -> gt_eval_v7/rediportal/scoring):
#   none       gt_eval_v6/prewindow/results_<tier>  *_step4_cov_edit_filtered
#   mismatch   gt_eval_v6/prewindow/results_<tier>  *_step5_nomm
#   clustered  gt_eval_v6/nomm/results_<tier>       *_pre100
# Output: gt_eval_v7/rediportal/scoring_<state>${OUT_SFX}/
#   for s in none mismatch clustered; do qsub human/pm/scripts/run_gt_rediportal_v7_states.sh $s; done
# v7 REDIportal truths (user, 2026-10-12): PGP1 SNV positions, dbSNP common and multiallelic sites
# removed, then the cov10_alt1 / alt3 / alt3_e05 gates. PGP1 removal is the default here
# (GT_REDI_EXCLUDE_BED; -v GT_REDI_EXCLUDE_BED= keeps PGP1 sites), multiallelic removal the
# default of gt_score_rediportal_v6.R (GT_REDI_DROP_MULTI=0 keeps them).
#$ -S /bin/bash
#$ -N gt_redi_v7s
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
set -euo pipefail
S=$1
B=${SCRNA_DIR:-.}
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
case $S in
  none)      TREE=prewindow; RE='_step4_cov_edit_filtered\.tsv\.gz$' ;;
  mismatch)  TREE=prewindow; RE='_step5_nomm\.tsv\.gz$' ;;
  clustered) TREE=nomm;      RE='_pre100\.tsv\.gz$' ;;
  *) echo "unknown state $S"; exit 1 ;;
esac
export GT_REDI_EXCLUDE_BED=${GT_REDI_EXCLUDE_BED-$B/human/pm/pgp1_snps/PGP1_hu43860C_hg38_snvs.bed}
OUT=$B/human/pm/gt_eval_v7/rediportal/scoring_$S${OUT_SFX:-}
mkdir -p "$OUT"
echo "OUT=$OUT GT_REDI_EXCLUDE_BED=$GT_REDI_EXCLUDE_BED GT_REDI_DROP_MULTI=${GT_REDI_DROP_MULTI:-1}"
GT_RES_TEMPLATE="gt_eval_v6/$TREE/results_%s" GT_CALL_RE="$RE" GT_REDI_OUT=$OUT \
  $RS human/scripts/gt_score_rediportal_v6.R
echo "done $S"
