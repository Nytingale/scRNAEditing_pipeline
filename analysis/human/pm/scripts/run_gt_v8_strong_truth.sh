#!/bin/bash
# run_gt_v8_strong_truth.sh   (2026-10-13)
# Scores the strong-NTC-evidence sensitivity truths (gt_v8_strong_truth.R) on the v8 call state (mismatch),
# with the v8 per-unit recall denominator and the 11-cell-type tier, then ranks (gt_rank_v7.R).
#   tasks 1-3 strong_alt10 x arm/sample/celltype, 4-6 strong_alt20
#   Rscript human/scripts/gt_v8_strong_truth.R 10 20
#   qsub -t 1-6 human/pm/scripts/run_gt_v8_strong_truth.sh
#   then: for a in 10 20; do Rscript human/scripts/gt_rank_v7.R gt_eval_v8/sensitivity/strong_alt$a/scoring/mismatch; done
#$ -S /bin/bash
#$ -N gt_v8_strong
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -l h_rt=4:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$B"
ALTS=(10 20)
TIERS=(arm sample celltype)
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
i=$(( ${SGE_TASK_ID:-1} - 1 ))
alt=${ALTS[$(( i / 3 ))]}; tier=${TIERS[$(( i % 3 ))]}
D=$V8/sensitivity/strong_alt$alt
echo "alt=$alt tier=$tier"
GT_V7_TRUTH=$D/truth/gt_v7_truth_sites.tsv.gz GT_UNIT_COV=$V8/unit_cov/ko_truth_unit_cov.tsv.gz \
GT_UNIT_KEEP=${GT_UNIT_KEEP:-$V8/unit_cov/celltype_keep_min20.txt} \
  $RS human/scripts/gt_score_v7.R "$tier" "$PM/${PUB[$tier]}" "$D/${SCN:-scoring}/mismatch/$tier" '_step5_nomm\.tsv\.gz$'
echo "done"
