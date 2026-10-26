#!/bin/bash
# mmonly_swap_v8.sh   (2026-10-15)
# After mmonly_rerun_v8.sh: every manifest row must have exit 0. Moves each published *_step5_nomm.tsv.gz (joint
# w1 x w2 window search) to gt_eval_v8/_superseded/step5_joint_w1w2_2026-10-15/<same relative path> and copies the
# w1-only step-5 table (gt_eval_v8/mmonly/results/...) in its place, so every scorer keeps reading the same paths.
# Also refreshes the replicate k-filter sets (replicates_kfilter_v8.R reads the per-sample step-5 tables).
#   bash human/pm/scripts/mmonly_swap_v8.sh   (or qsub, as mmonly_rescore_v8.sh does)
#$ -S /bin/bash
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_rt=4:00:00
set -euo pipefail
PM=${SCRNA_DIR:-.}/human/pm
P=$PM/gt_eval_v8/mmonly
SUP=$PM/gt_eval_v8/_superseded/step5_joint_w1w2_2026-10-15
bad=0; n=0
while IFS=$'\t' read -r tier wd pub cmd; do
  base=$(basename "$pub" _pre100.tsv.gz)
  st=$P/status/${tier}__${base}.tsv
  [ -s "$st" ] && [ "$(cut -f3 "$st")" = 0 ] || { echo "NOT DONE: $tier $base" >&2; bad=$((bad + 1)); }
done < <(tail -n +2 "$P/manifest.tsv")
# exit 100: under SGE the job goes to error state and every -hold_jid dependant stays held
[ "$bad" -eq 0 ] || { echo "ERROR: $bad rows not done; nothing swapped" >&2; exit 100; }

while IFS=$'\t' read -r tier wd pub cmd; do
  base=$(basename "$pub" _pre100.tsv.gz)
  rel=$(dirname "${pub#"$PM"/}")
  new=$P/results/$rel/${base}_step5_nomm.tsv.gz
  old=$PM/$rel/${base}_step5_nomm.tsv.gz
  [ -e "$new" ] || { echo "missing $new" >&2; exit 1; }
  if [ -e "$old" ] && [ ! -e "$SUP/$rel/${base}_step5_nomm.tsv.gz" ]; then
    mkdir -p "$SUP/$rel"; mv "$old" "$SUP/$rel/"
  fi
  cp "$new" "$old"
  n=$((n + 1))
done < <(tail -n +2 "$P/manifest.tsv")
echo "swapped $n step-5 tables (old copies: $SUP)"
rm -f "$PM/gt_eval_v8/replicates/results_k2_arm/.done"
