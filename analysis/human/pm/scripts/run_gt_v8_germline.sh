#!/bin/bash
# run_gt_v8_germline.sh  [submit|setup|collect|score]   (2026-10-13)
# Germline false-positive rate (user 2026-10-13): unit coverage at the germline sites (gt_v8_germline_sites.R)
# in the 28 NTC unit BAMs (gt_v8_unit_cov.sh, OUT = gt_eval_v8/germline/unit_cov), then
# gt_score_germline_v8.R per tier on the mismatch state (task 1-3 = arm, sample, celltype).
#   HOLD=<gl_sites job> bash run_gt_v8_germline.sh submit
#$ -S /bin/bash
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
G=$PM/gt_eval_v8/germline
SELF=$B/human/pm/scripts/run_gt_v8_germline.sh
UC=$B/human/pm/scripts/gt_v8_unit_cov.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export OUT=$G/unit_cov TRUTH=$G/germline_sites.tsv.gz
STAGE=${STAGE:-${1:-submit}}
cd "$B"

if [ "$STAGE" = "submit" ]; then
  V="OUT=$OUT,TRUTH=$TRUTH"
  j1=$(qsub -terse -N gl_setup ${HOLD:+-hold_jid $HOLD} -v STAGE=setup "$SELF")
  j2=$(qsub -terse -N gl_ucov -t 1-28 -hold_jid "$j1" -v "$V" "$UC" | cut -d. -f1)
  j3=$(qsub -terse -N gl_collect -hold_jid "$j2" -v STAGE=collect "$SELF")
  j4=$(qsub -terse -N gl_score -t 1-3 -l h_vmem=48G -hold_jid "$j3" -v STAGE=score "$SELF" | cut -d. -f1)
  echo "setup $j1 | unit cov $j2 (1-28) | collect $j3 | score $j4 (1-3)"
  exit 0
fi
if [ "$STAGE" = "setup" ]; then bash "$UC" setup; [ "$(wc -l < "$OUT/tasks.tsv")" -eq 28 ]; exit 0; fi
if [ "$STAGE" = "collect" ]; then bash "$UC" collect; exit 0; fi
if [ "$STAGE" = "score" ]; then
  TIERS=(arm sample celltype); tier=${TIERS[$(( SGE_TASK_ID - 1 ))]}
  GT_UNIT_KEEP=${GT_UNIT_KEEP:-$PM/gt_eval_v8/unit_cov/celltype_keep_min20.txt} $RS human/scripts/gt_score_germline_v8.R "$tier"
  exit 0
fi
echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
