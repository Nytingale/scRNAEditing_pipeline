#!/bin/bash
# run_gt_v8_rescore.sh   (2026-10-13)
# v8 rescoring with per-unit recall denominators (user 2026-10-13), written into gt_eval_v8/.
#   tasks 1-12   KO truth (gt_score_v7.R, GT_UNIT_COV = gt_eval_v8/unit_cov/ko_truth_unit_cov.tsv.gz),
#                state x tier -> gt_eval_v8/scoring/<state>/<tier>/
#   tasks 13-16  REDIportal (gt_score_rediportal_v8.R: per-truth universe, per-unit denominator),
#                state -> gt_eval_v8/rediportal/scoring_<state>/
# Call-set window states as run_gt_v7_variants_score.sh / run_gt_rediportal_v7_states.sh.
#   tasks 17-19  KO truth with the region mask (gt_eval_v7 variant v7_mask_nomm_pgp1, a subset of the v8 positions),
#                mismatch state x tier -> gt_eval_v8/mask_compare/mismatch/<tier>/ (GTDEF_v8_mask_vs_nomm_plot.R panel d)
#   qsub -t 1-19 human/pm/scripts/run_gt_v8_rescore.sh
# Then: gt_rank_v7.R per state and gt_v8_compare.R.
#$ -S /bin/bash
#$ -N gt_v8_rescore
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=48G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
# cell-type tier: the cell types with >= 20 truth sites (>= 10 reads) under all nine truths (user 2026-10-13;
# unit_cov/celltype_covered_positives.tsv, TERATOMA_v8_density_celltype_truth_sites); -v GT_UNIT_KEEP= uses all 23
export GT_UNIT_KEEP=${GT_UNIT_KEEP-$V8/unit_cov/celltype_keep_min20.txt}
# -v GT_UNIT_RULE=gate SC=$V8/scoring_pw RD=$V8/rediportal_pw: power-gated positives need unit coverage x merged NTC
# editing level >= 3 (gt_score_v7.R); depth-gated and REDIportal keep >= 10 reads
export GT_UNIT_RULE=${GT_UNIT_RULE:-cov}
SC=${SC:-$V8/scoring}
RD=${RD:-$V8/rediportal}
echo "GT_UNIT_KEEP=$GT_UNIT_KEEP GT_UNIT_RULE=$GT_UNIT_RULE SC=$SC RD=$RD"
STATES=(none mismatch clustered both)
TIERS=(arm sample celltype)
# Call sets (2026-10-15): the v8 caller runs (run_pipeline_pm_gt_v8.sh; step-4 / step-5 tables next to each pre100)
# and their no-window replay (nomm_rerun_v8.sh -> *_v8_nomm). Before: v6 runs + gt_eval_v6/prewindow, v6_nomm.
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
RES_DIRS="arm=${PUB[arm]},sample=${PUB[sample]},celltype=${PUB[celltype]}"
i=$(( ${SGE_TASK_ID:-1} - 1 ))
if [ $i -ge 16 ]; then
  tier=${TIERS[$(( i - 16 ))]}
  echo "KO v7_mask_nomm_pgp1 state=mismatch tier=$tier"
  GT_V7_TRUTH=$PM/gt_eval_v7/variants/truth/v7_mask_nomm_pgp1/gt_v7_truth_sites.tsv.gz \
  GT_UNIT_COV=$V8/unit_cov/ko_truth_unit_cov.tsv.gz \
    $RS human/scripts/gt_score_v7.R "$tier" "$PM/${PUB[$tier]}" "${MC:-$V8/mask_compare}/mismatch/$tier" \
    '_step5_nomm\.tsv\.gz$'
elif [ $i -lt 12 ]; then
  state=${STATES[$(( i / 3 ))]}; tier=${TIERS[$(( i % 3 ))]}
  case $state in
    none)      res=$PM/${PUB[$tier]};      re='_step4_cov_edit_filtered\.tsv\.gz$' ;;
    mismatch)  res=$PM/${PUB[$tier]};      re='_step5_nomm\.tsv\.gz$' ;;
    clustered) res=$PM/${PUB[$tier]}_nomm; re='_pre100\.tsv\.gz$' ;;
    both)      res=$PM/${PUB[$tier]};      re='_pre100\.tsv\.gz$' ;;
  esac
  echo "KO state=$state tier=$tier res=$res"
  GT_V7_TRUTH=$V8/truth/gt_v7_truth_sites.tsv.gz GT_UNIT_COV=$V8/unit_cov/ko_truth_unit_cov.tsv.gz \
    $RS human/scripts/gt_score_v7.R "$tier" "$res" "$SC/$state/$tier" "$re"
else
  state=${STATES[$(( i - 12 ))]}
  dirs=$RES_DIRS
  case $state in
    none)      RE='_step4_cov_edit_filtered\.tsv\.gz$' ;;
    mismatch)  RE='_step5_nomm\.tsv\.gz$' ;;
    clustered) RE='_pre100\.tsv\.gz$'; dirs="arm=${PUB[arm]}_nomm,sample=${PUB[sample]}_nomm,celltype=${PUB[celltype]}_nomm" ;;
    both)      RE='_pre100\.tsv\.gz$' ;;
  esac
  echo "REDIportal state=$state dirs=$dirs"
  GT_REDI_EXCLUDE_BED=$PM/pgp1_snps/PGP1_hu43860C_hg38_snvs.bed GT_RES_DIRS="$dirs" GT_CALL_RE="$RE" \
    GT_REDI_OUT=$RD/scoring_$state $RS human/scripts/gt_score_rediportal_v8.R
fi
echo "done"
