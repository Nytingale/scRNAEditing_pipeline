#!/bin/bash
# =============================================================================
# run_gt_v8_ko_floor.sh  [submit|eval|truth|ucov_setup|ucov|ucov_collect|score|rank]   (2026-10-17)
#
# KO coverage floor sweep (user 2026-10-17): v8 truth sets built exactly as the depth gate (NTC cov >= 10 and
# KO cov >= 10 per arm) with the KO floor at 5 (gate depth5) and 3 (depth3); NTC rules unchanged (rule 2
# already implies >= 20 merged NTC reads). Same chain as the v8 truth, reusing the v6 spectrum12 pileups:
#   eval        gt_spectrum12_eval_v6.R  GT_COVER_GATE=depth<N> -> gt_eval_v6/spectrum12/eval_depth<N>/
#   truth       gt_truth_agtc_test_v6.R (GT_TRUTH_ONLY, GT_GATES) -> gt_truth_v7.R -> gt_truth_v7_variants.R
#               (nearby-mismatch w1 searched per gate, as the v8 truth) -> ko_floor/truth/v7_nomm_pgp1/
#   ucov        gt_v8_unit_cov.sh on the new truth positions -> ko_floor/unit_cov/ko_truth_unit_cov.tsv.gz
#   score       gt_score_v7.R, mismatch state (*_step5_nomm), per tier; cell-type tier = celltype_keep_min20 (as v8)
#   rank        gt_rank_v7.R -> ko_floor/scoring/mismatch/config_rank.tsv; gt_v8_ko_floor_compare.R
# Outputs: human/pm/gt_eval_v8/ko_floor/
#   bash human/pm/scripts/run_gt_v8_ko_floor.sh submit
#$ -S /bin/bash
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
KF=$V8/ko_floor
SELF=$PM/scripts/run_gt_v8_ko_floor.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:${PATH}"
FLOORS=(5 3)
GATE_LIST=$(printf "depth%s," "${FLOORS[@]}"); GATE_LIST=${GATE_LIST%,}
export GT_SP12_DIR=$PM/gt_eval_v6/spectrum12 GT_AGTC_DIR=$KF/agtc_test GT_V7_TRUTH_DIR=$KF/truth_v7 \
       GT_V7_VARIANTS_DIR=$KF/truth GT_GATES=$GATE_LIST GT_TRUTH_ONLY=1
TRUTH=$KF/truth/v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
TIERS=(arm sample celltype)
mkdir -p "$KF" "$B/logs"
cd "$B"
STAGE="${1:-submit}"

if [ "$STAGE" = "submit" ]; then
  q="qsub -terse -q slice"
  je=$($q -N kf_eval -l h_vmem=96G -t 1-${#FLOORS[@]} "$SELF" eval | cut -d. -f1)
  jt=$($q -N kf_truth -l h_vmem=96G -hold_jid "$je" "$SELF" truth)
  js=$($q -N kf_ucov_setup -hold_jid "$jt" "$SELF" ucov_setup)
  nb=$(ls $PM/02_gt_bams_umi_v6/merged_NTC_ctrl.bam $PM/02_gt_bams_umi_v6/T[1-4]_NTC_ctrl.bam \
          $PM/02_gt_bams_celltype_umi_v6/NTC_ctrl__*.bam | wc -l)
  ju=$($q -N kf_ucov -l h_vmem=4G -t 1-"$nb" -hold_jid "$js" "$SELF" ucov | cut -d. -f1)
  jc=$($q -N kf_ucov_collect -hold_jid "$ju" "$SELF" ucov_collect)
  jsc=$($q -N kf_score -l h_vmem=48G -l h_rt=12:00:00 -t 1-3 -hold_jid "$jc" "$SELF" score | cut -d. -f1)
  jr=$($q -N kf_rank -l h_vmem=16G -hold_jid "$jsc" "$SELF" rank)
  echo "eval $je | truth $jt | ucov_setup $js | ucov $ju (1-$nb) | collect $jc | score $jsc | rank $jr"
  exit 0
fi

if [ "$STAGE" = "eval" ]; then
  f=${FLOORS[$(( SGE_TASK_ID - 1 ))]}
  GT_COVER_GATE=depth$f $RS human/scripts/gt_spectrum12_eval_v6.R
  exit 0
fi

if [ "$STAGE" = "truth" ]; then
  for f in "${FLOORS[@]}"; do [ -s "$GT_SP12_DIR/eval_depth$f/truth12_sites.tsv.gz" ] || { echo "ERROR: eval_depth$f missing" >&2; exit 1; }; done
  $RS human/scripts/gt_truth_agtc_test_v6.R
  $RS human/scripts/gt_truth_v7.R
  $RS human/scripts/gt_truth_v7_variants.R
  exit 0
fi

UC_ENV="OUT=$KF/unit_cov TRUTH=$TRUTH"
if [ "$STAGE" = "ucov_setup" ]; then env $UC_ENV bash "$PM/scripts/gt_v8_unit_cov.sh" setup; exit 0; fi
if [ "$STAGE" = "ucov" ]; then env $UC_ENV bash "$PM/scripts/gt_v8_unit_cov.sh"; exit 0; fi
if [ "$STAGE" = "ucov_collect" ]; then env $UC_ENV bash "$PM/scripts/gt_v8_unit_cov.sh" collect; exit 0; fi

if [ "$STAGE" = "score" ]; then
  tier=${TIERS[$(( SGE_TASK_ID - 1 ))]}
  GT_V7_TRUTH=$TRUTH GT_UNIT_COV=$KF/unit_cov/ko_truth_unit_cov.tsv.gz GT_UNIT_KEEP=${GT_UNIT_KEEP:-$V8/unit_cov/celltype_keep_min20.txt} \
    $RS human/scripts/gt_score_v7.R "$tier" "$PM/${PUB[$tier]}" "$KF/${SCN:-scoring}/mismatch/$tier" '_step5_nomm\.tsv\.gz$'
  exit 0
fi

if [ "$STAGE" = "rank" ]; then
  $RS human/scripts/gt_rank_v7.R "$KF/${SCN:-scoring}/mismatch"
  $RS human/scripts/gt_v8_ko_floor_compare.R
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
