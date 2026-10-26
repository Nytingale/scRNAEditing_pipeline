#!/bin/bash
# run_gt_split_score.sh   (2026-10-13)
# Split-by-cells test, scoring. Half-B calls (results_gt_split_B_v8, mismatch state = *_step5_nomm, the main
# v8 state) scored per tier with gt_score_v7.R (per-unit recall: positives with >= 10 reads in the half-B unit;
# cell-type tier = the 11 kept cell types):
#   tasks 1-3  against the half-A truth  (split_cells/truth/v7_nomm_pgp1; reads disjoint from the calls)
#              -> split_cells/scoring/truthA/<tier>/
#   tasks 4-6  against the full v8 truth (all NTC cells, so it contains the half-B reads)
#              -> split_cells/scoring/truthFull/<tier>/
# Unit coverage: gt_v8_unit_cov.sh on the half-B UMI BAMs at the union of both truths' positions:
#   OUT=$S/unit_cov TRUTH="<A> <full>" BAMS="$S/bams/umi/*_B.bam" bash gt_v8_unit_cov.sh setup
#   qsub -t 1-16 -v OUT=...,TRUTH=...,BAMS=... gt_v8_unit_cov.sh ; ... collect
#   qsub -t 1-6 run_gt_split_score.sh ; then Rscript human/scripts/gt_split_compare.R
#$ -S /bin/bash
#$ -N gt_split_score
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=48G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
S=$PM/gt_eval_v8/split_cells
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
export GT_UNIT_KEEP=${GT_UNIT_KEEP:-$PM/gt_eval_v8/unit_cov/celltype_keep_min20.txt} GT_UNIT_COV=$S/unit_cov/ko_truth_unit_cov.tsv.gz
TIERS=(arm sample celltype)
i=$(( ${SGE_TASK_ID:-1} - 1 ))
tier=${TIERS[$(( i % 3 ))]}
if [ $i -lt 3 ]; then tn=truthA; tr=$S/truth/v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz
else tn=truthFull; tr=$PM/gt_eval_v8/truth/gt_v7_truth_sites.tsv.gz; fi
echo "half-B calls, $tn, tier $tier"
GT_V7_TRUTH=$tr $RS human/scripts/gt_score_v7.R "$tier" "$PM/results_gt_split_B_v8" "$S/${SCN:-scoring}/$tn/$tier" '_step5_nomm\.tsv\.gz$'
echo done
