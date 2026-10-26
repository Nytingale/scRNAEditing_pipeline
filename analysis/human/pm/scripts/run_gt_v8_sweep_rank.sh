#!/bin/bash
# run_gt_v8_sweep_rank.sh  [submit|truth|ucov_setup|ucov_collect|score|rank]   (2026-10-14)
# Rank sensitivity to every v8 truth threshold (user 2026-10-14): one threshold at a time, +-1 step around
# the centre (gt_robustness_v8.R, GT_ROBUST_TRUTH_OUT; depletion excluded: it is the t25 / t50 / t75 axis),
# each neighbour truth scored like the main v8 ranking (call state mismatch, per-unit recall, 11 cell types)
# and ranked (gt_rank_v7.R). The centre is rebuilt and scored too (must reproduce the main ranking).
#   truth         -> gt_eval_v8/sweep_rank/truths/<axis>_<value>/gt_v7_truth_sites.tsv.gz, truth_counts.tsv
#   ucov_*        unit coverage at the union of all neighbour truth positions -> sweep_rank/unit_cov/
#   score         array 1..(3 x truths): truth x tier -> sweep_rank/scoring/<truth>/mismatch/<tier>/
#   rank          gt_rank_v7.R per truth -> .../mismatch/config_rank.tsv
#   bash run_gt_v8_sweep_rank.sh submit     (chains everything with -hold_jid)
#   then: Rscript human/scripts/gt_v8_sweep_rank_compare.R
#$ -S /bin/bash
#$ -N gt_sweep_rank
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
W=$V8/sweep_rank
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
SELF=$PM/scripts/run_gt_v8_sweep_rank.sh
UCOV=$PM/scripts/gt_v8_unit_cov.sh
TIERS=(arm sample celltype)
# v8 caller runs (2026-10-14; BCFtools BAQ on, REF=N fixed); the v6 call sets were gt_eval_v6/prewindow/results_<tier>
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
N_TRUTH=15   # centre + 14 neighbours (7 axes x 2)
cd "$B"
STAGE=${1:-${STAGE:-}}
SCN=${SCN:-scoring}

if [ "$STAGE" = "submit" ]; then
  jt=$(qsub -terse -N swr_truth "$SELF" truth)
  js=$(qsub -terse -N swr_usetup -l h_vmem=8G -hold_jid "$jt" "$SELF" ucov_setup)
  ju=$(qsub -terse -N swr_ucov -t 1-28 -l h_vmem=4G -hold_jid "$js" -v OUT="$W/unit_cov" "$UCOV" | cut -d. -f1)
  jc=$(qsub -terse -N swr_ucollect -l h_vmem=8G -hold_jid "$ju" "$SELF" ucov_collect)
  jsc=$(qsub -terse -N swr_score -t 1-$(( 3 * N_TRUTH )) -hold_jid "$jc" "$SELF" score | cut -d. -f1)
  jr=$(qsub -terse -N swr_rank -l h_vmem=8G -hold_jid "$jsc" "$SELF" rank)
  echo "truth $jt | ucov setup $js | ucov $ju (1-28) | collect $jc | score $jsc (1-$(( 3 * N_TRUTH ))) | rank $jr"
  exit 0
fi

if [ "$STAGE" = "truth" ]; then
  rm -rf "$W/truths"; mkdir -p "$W/truths"
  GT_ROBUST_TRUTH_OUT=$W/truths /usr/bin/time -v $RS human/scripts/gt_robustness_v8.R
  n=$(ls -d "$W"/truths/*/ | wc -l)
  [ "$n" -eq "$N_TRUTH" ] || { echo "ERROR: $n truths, expected $N_TRUTH" >&2; exit 1; }
  exit 0
fi

if [ "$STAGE" = "ucov_setup" ]; then
  OUT=$W/unit_cov TRUTH="$(ls "$W"/truths/*/gt_v7_truth_sites.tsv.gz | tr '\n' ' ')" bash "$UCOV" setup
  [ "$(wc -l < "$W/unit_cov/tasks.tsv")" -eq 28 ] || { echo "ERROR: expected 28 unit BAMs" >&2; exit 1; }
  exit 0
fi

if [ "$STAGE" = "ucov_collect" ]; then
  n=$(ls "$W"/unit_cov/ko_truth/*.tsv.gz | wc -l)
  [ "$n" -eq 28 ] || { echo "ERROR: $n / 28 unit coverage files" >&2; exit 1; }
  OUT=$W/unit_cov bash "$UCOV" collect
  exit 0
fi

if [ "$STAGE" = "score" ]; then
  i=$(( SGE_TASK_ID - 1 ))
  t=$(ls -d "$W"/truths/*/ | sort | sed -n "$(( i / 3 + 1 ))p"); t=$(basename "$t")
  tier=${TIERS[$(( i % 3 ))]}
  echo "truth $t tier $tier"
  GT_V7_TRUTH=$W/truths/$t/gt_v7_truth_sites.tsv.gz GT_UNIT_COV=$W/unit_cov/ko_truth_unit_cov.tsv.gz \
  GT_UNIT_KEEP=${GT_UNIT_KEEP:-$V8/unit_cov/celltype_keep_min20.txt} \
    $RS human/scripts/gt_score_v7.R "$tier" "$PM/${PUB[$tier]}" "$W/$SCN/$t/mismatch/$tier" '_step5_nomm\.tsv\.gz$'
  echo done
  exit 0
fi

if [ "$STAGE" = "rank" ]; then
  for d in "$W"/$SCN/*/mismatch; do
    for tier in "${TIERS[@]}"; do [ -s "$d/$tier/gt_v7_scores_reported.tsv" ] || { echo "ERROR: missing $d/$tier" >&2; exit 1; }; done
    $RS human/scripts/gt_rank_v7.R "$d"
  done
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
