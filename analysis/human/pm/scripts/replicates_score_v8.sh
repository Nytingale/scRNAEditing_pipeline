#!/bin/bash
# replicates_score_v8.sh   (2026-10-15; Reviewer 3 replicates; held on replicates_jacusa2_v8.sh; from _v6)
# Arm-level NTC call sets (main state, _step5_nomm) under the replicate strategies, scored against the three
# REDIportal truths (gt_score_rediportal_v8.R) and the six KO truths (gt_score_v7.R, arm tier):
#   pooled       results_gt_v8 (merged NTC arm; existing v8 scoring, rediportal/scoring_mismatch, scoring/mismatch/arm)
#   k2 / k3 / k4 in >= k of the T1-T4 per-sample call sets (replicates_kfilter_v8.R), every caller
#   jacusa2rep   JACUSA2 call-1 with T1-T4 as replicates (replicates_jacusa2_v8.sh)
# Then replicates_summary_v8.R.
# Output: gt_eval_v8/replicates/scoring/{rediportal,ko}_<strategy>/, gt_eval_v8/replicates/summary/
#$ -S /bin/bash
#$ -N rep_score8
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -pe threaded 2
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
RR=$V8/replicates
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
unset GT_UNIT_KEEP
[ -s "$RR/results_k2_arm/.done" ] || { $RS human/scripts/replicates_kfilter_v8.R && touch "$RR/results_k2_arm/.done"; }
mkdir -p "$RR/_empty"
for s in ${STRATS:-k2 k3 k4 jacusa2rep}; do
  res=gt_eval_v8/replicates/results_${s}_arm
  GT_REDI_EXCLUDE_BED=$PM/pgp1_snps/PGP1_hu43860C_hg38_snvs.bed \
  GT_RES_DIRS="arm=$res,sample=gt_eval_v8/replicates/_empty,celltype=gt_eval_v8/replicates/_empty" \
  GT_CALL_RE='_step5_nomm\.tsv\.gz$' GT_REDI_OUT=$RR/scoring/rediportal_$s \
    $RS human/scripts/gt_score_rediportal_v8.R
  GT_V7_TRUTH=$V8/truth/gt_v7_truth_sites.tsv.gz GT_UNIT_COV=$V8/unit_cov/ko_truth_unit_cov.tsv.gz \
    $RS human/scripts/gt_score_v7.R arm "$PM/$res" "$RR/scoring/ko_$s" '_step5_nomm\.tsv\.gz$'
done
[ -s "$RR/scoring/rediportal_jacusa2rep/unit_scores.tsv" ] && $RS human/scripts/replicates_summary_v8.R
echo done
