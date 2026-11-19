#!/bin/bash
# run_gt_v9_germline.sh  [submit|sites|setup|collect|score]   (2026-11-17)
# v9 copy of run_gt_v8_germline.sh (R3.9 homozygous PGP1 variants; response paragraph 496): germline sites rebuilt
# from the seed-8787 spectrum12 pileups (gt_v8_germline_sites.R, disjoint from the v9 KO positives), unit coverage
# in the 28 seeded NTC unit BAMs (gt_v8_unit_cov.sh), then gt_score_germline_v8.R per tier for the 76 ranked
# configurations: published trees (58), Cellsnp-lite hybrid (gt_eval_v8/cellsnp_hybrid/results/centre, 2) and the
# UMI-aware MarkDuplicates trees (16). Output: gt_eval_v9/germline/scoring_pwcov/<set>/<tier>/.
#   bash human/pm/scripts/run_gt_v9_germline.sh submit   (on the Slurm login node, cluster/slurm on PATH)
#$ -S /bin/bash
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
export GT_GERM_DIR=$PM/gt_eval_v9/germline
G=$GT_GERM_DIR
SELF=$PM/scripts/run_gt_v9_germline.sh
UC=$PM/scripts/gt_v8_unit_cov.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export OUT=$G/unit_cov TRUTH=$G/germline_sites.tsv.gz
STAGE=${STAGE:-${1:-submit}}
cd "$B"

case $STAGE in
submit)
  mkdir -p "$G"
  j0=$(qsub -terse -N gl9_sites -l h_vmem=32G -v STAGE=sites "$SELF")
  j1=$(qsub -terse -N gl9_setup -hold_jid "$j0" -v STAGE=setup "$SELF")
  j2=$(qsub -terse -N gl9_ucov -t 1-28 -hold_jid "$j1" -v "OUT=$OUT,TRUTH=$TRUTH" "$UC" | cut -d. -f1)
  j3=$(qsub -terse -N gl9_collect -hold_jid "$j2" -v STAGE=collect "$SELF")
  j4=$(qsub -terse -N gl9_score -t 1-9 -l h_vmem=48G -hold_jid "$j3" -v STAGE=score "$SELF" | cut -d. -f1)
  echo "sites $j0 | setup $j1 | unit cov $j2 (1-28) | collect $j3 | score $j4 (1-9)"
  ;;
sites)
  GT_GERM_TRUTH=$PM/gt_eval_v9/truth/gt_v9_positive_sites.tsv.gz $RS human/scripts/gt_v8_germline_sites.R
  ;;
setup) bash "$UC" setup; [ "$(wc -l < "$OUT/tasks.tsv")" -eq 28 ] ;;
collect) bash "$UC" collect ;;
score)
  TIERS=(arm sample celltype); i=$(( SGE_TASK_ID - 1 )); tier=${TIERS[$(( i % 3 ))]}
  declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
  case $(( i / 3 )) in
    0) set=pub; res=$PM/${PUB[$tier]} ;;
    1) set=ch; res=$PM/gt_eval_v8/cellsnp_hybrid/results/centre/${PUB[$tier]} ;;
    2) set=umimd; res=$PM/${PUB[$tier]}_umimd ;;
  esac
  GT_UNIT_KEEP=$PM/gt_eval_v8/unit_cov/celltype_keep_pwcov.txt \
    $RS human/scripts/gt_score_germline_v8.R "$tier" "$res" "$G/scoring_pwcov/$set/$tier"
  ;;
*) echo "ERROR: unknown stage '$STAGE'" >&2; exit 1 ;;
esac
