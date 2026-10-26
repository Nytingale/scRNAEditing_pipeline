#!/bin/bash
# run_gt_v7_variants_score.sh   (2026-10-10)
# One array task = one call-set window state x tier; scores the v7 truth variants (v7, v7_mask, v7_nomm, v7_mask_nomm,
# each also _pgp1; subset with -v VARIANTS="v7_pgp1 ...")
# (gt_eval_v7/variants/truth/<variant>/) with gt_score_v7.R.
# Window states (TERATOMA_window_state_ranksum_plot.R):
#   none       no window filter          gt_eval_v6/prewindow/results_<tier>  *_step4_cov_edit_filtered
#   mismatch   nearby-mismatch pass only gt_eval_v6/prewindow/results_<tier>  *_step5_nomm
#   clustered  clustered-sites pass only results_*_v6_nomm                    *_pre100
#   both       both passes (published)   results_*_v6                         *_pre100
# Output: gt_eval_v7/variants/scoring/<state>/<variant>/<tier>/ (or $SC_OUT/<state>/<variant>/<tier>/)
#   qsub -t 1-12 [-v SC_OUT=<dir>] human/pm/scripts/run_gt_v7_variants_score.sh
#$ -S /bin/bash
#$ -N gt_v7var_sc
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -l h_rt=8:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
SC_OUT=${SC_OUT:-$PM/gt_eval_v7/variants/scoring}
cd "$B"
STATES=(none mismatch clustered both)
TIERS=(arm sample celltype)
i=$(( ${SGE_TASK_ID:-1} - 1 ))
state=${STATES[$(( i / 3 ))]}
tier=${TIERS[$(( i % 3 ))]}
declare -A PUB=([arm]=results_gt_v6 [sample]=results_gt_persample_v6 [celltype]=results_gt_celltype_v6)
case $state in
  none)      res=$PM/gt_eval_v6/prewindow/results_$tier; re='_step4_cov_edit_filtered\.tsv\.gz$' ;;
  mismatch)  res=$PM/gt_eval_v6/prewindow/results_$tier; re='_step5_nomm\.tsv\.gz$' ;;
  clustered) res=$PM/${PUB[$tier]}_nomm;                 re='_pre100\.tsv\.gz$' ;;
  both)      res=$PM/${PUB[$tier]};                      re='_pre100\.tsv\.gz$' ;;
esac
echo "state=$state tier=$tier res=$res"
for v in ${VARIANTS:-v7 v7_mask v7_nomm v7_mask_nomm v7_pgp1 v7_mask_pgp1 v7_nomm_pgp1 v7_mask_nomm_pgp1}; do
  GT_V7_TRUTH=$PM/gt_eval_v7/variants/truth/$v/gt_v7_truth_sites.tsv.gz \
    $RS human/scripts/gt_score_v7.R "$tier" "$res" "$SC_OUT/$state/$v/$tier" "$re"
done
echo "done"
