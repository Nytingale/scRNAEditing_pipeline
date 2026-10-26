#!/bin/bash
# replicates_reachable_v8.sh   (2026-10-15; Reviewer 3 replicates, fairness check)
# Pooled vs replicate strategies with recall restricted to positives with >= 10 reads in each of T1-T4 NTC
# (human/scripts/replicates_reachable_v8.R). KO truths rescored with gt_score_v7.R on the reachable coverage table.
# qsub -v REACH_MIN_TER=2 ... : >= 10 reads in at least 2 of T1-T4 instead (output reachable_min2/).
# Output: gt_eval_v8/replicates/reachable/
#$ -S /bin/bash
#$ -N rep_reach8
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -pe threaded 2
#$ -l h_rt=6:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
RR=$V8/replicates
export REACH_MIN_TER=${REACH_MIN_TER:-4}
if [ "$REACH_MIN_TER" = 4 ]; then O=$RR/reachable; else O=$RR/reachable_min$REACH_MIN_TER; fi
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
unset GT_UNIT_KEEP
$RS human/scripts/replicates_reachable_v8.R prep
for s in pooled k2 k3 k4 jacusa2rep; do
  if [ "$s" = pooled ]; then res=$PM/results_gt_v8; else res=$RR/results_${s}_arm; fi
  GT_V7_TRUTH=$V8/truth/gt_v7_truth_sites.tsv.gz GT_UNIT_COV=$O/ko_truth_unit_cov_reachable.tsv.gz \
    $RS human/scripts/gt_score_v7.R arm "$res" "$O/ko_$s" '_step5_nomm\.tsv\.gz$'
done
$RS human/scripts/replicates_reachable_v8.R score
echo done
