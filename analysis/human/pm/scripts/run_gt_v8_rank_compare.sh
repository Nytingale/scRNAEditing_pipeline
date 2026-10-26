#!/bin/bash
# run_gt_v8_rank_compare.sh   (2026-10-15)
# After run_gt_v8_rescore.sh (all 19 tasks): gt_rank_v7.R on gt_eval_v8/scoring/<state> and mask_compare/mismatch,
# then gt_v8_compare.R -> gt_eval_v8/compare_config_rank.tsv.
#$ -S /bin/bash
#$ -N gt_v8_rankcmp
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=4:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
V8=$B/human/pm/gt_eval_v8
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
for d in "$V8"/scoring/{none,mismatch,clustered,both} "$V8/mask_compare/mismatch"; do
  for t in arm sample celltype; do [ -s "$d/$t/gt_v7_scores.tsv" ] || { echo "ERROR: missing $d/$t" >&2; exit 1; }; done
  $RS human/scripts/gt_rank_v7.R "$d" > "$d/rank.log" 2>&1 || { tail "$d/rank.log" >&2; exit 1; }
  echo "ranked $d"
done
$RS human/scripts/gt_v8_compare.R
