#!/bin/bash
# aorta_umimd_score.sh   (user 2026-11-09; production run 2026-11-11)
# v9 scoring of the aorta UMI-aware MarkDuplicates call sets, same settings as aorta_v9_build.sh score (v9 truth and
# unit coverage, gate_cov, edit unit CTRL, GT_FP_RULE=ag, step5_nomm call tables).
# Default SRC=prod: results_v8_umimd (production, 16 configurations, bam type sampleUmiMd); array 1 arm, 2 sample.
# SRC=test: the old test sets (human/pm/scripts/umimd_run.sh TIER=aorta|aorta_ctrl, sampleUmiDup slot);
# array 1 arm, 2 sample (results_umimd_aorta), 3 arm_ctrl (results_umimd_aorta_ctrl: published UMI BAMs re-run,
# reproducibility control, scored as arm). The test sets have no step5_nomm tables (2026-11-11).
# Output: gt_eval_v9/umimd/<set>/gt_v7_scores.tsv ... ; then Rscript tmp/umimd_compare/umimd_compare.R aorta
# Usage: qsub -t 1-2 mouse/aorta_gse255661/scripts/aorta_umimd_score.sh   (SRC=test: qsub -v SRC=test -t 1-3 ...)
#$ -S /bin/bash
#$ -N aorta_umimd_score
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/gt_eval_v9/logs/
#$ -j y
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
A=$B/mouse/aorta_gse255661
V9=$A/gt_eval_v9
TR=$A/truth/eval_v9/excluded
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
export GT_FP_RULE=ag
SRC=${SRC:-prod}
SETS=(arm sample arm_ctrl)
set_name=${SETS[$(( SGE_TASK_ID - 1 ))]}
tier=${set_name%_ctrl}
if [ "$SRC" = prod ]; then
  [ "$set_name" = arm_ctrl ] && { echo "no arm_ctrl set for SRC=prod" >&2; exit 1; }
  res=$A/results_v8_umimd
else
  res=$( [ "$set_name" = arm_ctrl ] && echo "$A/results_umimd_aorta_ctrl" || echo "$A/results_umimd_aorta" )
fi
ure=$( [ "$tier" = arm ] && echo '^CTRL$' || echo '^Control_' )
cd "$B"
GT_V7_TRUTH=$TR/gt_aorta_truth_sites.tsv.gz GT_UNIT_COV=$TR/unit_cov.tsv.gz GT_UNIT_RULE=gate_cov \
GT_UNIT_EDIT_UNIT=CTRL GT_UNIT_RE="$ure" \
  Rscript human/scripts/gt_score_v7.R "$tier" "$res" "$V9/umimd/$set_name" '_step5_nomm\.tsv\.gz$'
