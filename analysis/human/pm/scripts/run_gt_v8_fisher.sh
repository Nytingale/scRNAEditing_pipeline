#!/bin/bash
# =============================================================================
# run_gt_v8_fisher.sh  [submit|eval|truth|fisher|ucov_setup|ucov|ucov_collect|score|rank]   (2026-10-18)
#
# Set B of HANDOFF_2026-10-17 4c-DECISION: v8 truth without the m = 3 gate, positives by a one-sided Fisher's
# exact test (KO vs merged NTC reads) plus the drop >= t. Same chain as run_gt_v8_ko_floor.sh:
#   eval        gt_spectrum12_eval_v6.R  GT_COVER_GATE=depth1 (NTC cov >= 10, KO cov >= 1) -> spectrum12/eval_depth1/
#   truth       gt_truth_agtc_test_v6.R -> gt_truth_v7.R -> gt_truth_v7_variants.R (w1 searched, as v8)
#               -> fisher/truth/v7_nomm_pgp1/
#   fisher      gt_v8_fisher_truth.R -> fisher/truth_fisher/gt_v7_truth_sites.tsv.gz, gates depth1 (effect size
#               only), fisher_bonf (arm p <= 0.025), fisher_bh (BH q <= 0.10); decoy / SNP chance positives
#   ucov        gt_v8_unit_cov.sh on the truth positions -> fisher/unit_cov/ko_truth_unit_cov.tsv.gz
#   score       gt_score_v7.R, mismatch state (*_step5_nomm), per tier; cell-type tier = celltype_keep_min20 (as v8)
#   rank        gt_rank_v7.R -> fisher/scoring/mismatch/config_rank.tsv; gt_v8_fisher_compare.R
# Outputs: human/pm/gt_eval_v8/fisher/
#   bash human/pm/scripts/run_gt_v8_fisher.sh submit
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
FT=$V8/fisher
SELF=$PM/scripts/run_gt_v8_fisher.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:${PATH}"
GATE=depth1
export GT_SP12_DIR=$PM/gt_eval_v6/spectrum12 GT_AGTC_DIR=$FT/agtc_test GT_V7_TRUTH_DIR=$FT/truth_v7 \
       GT_V7_VARIANTS_DIR=$FT/truth GT_GATES=$GATE GT_TRUTH_ONLY=1
TRUTH=$FT/truth_fisher/gt_v7_truth_sites.tsv.gz
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
TIERS=(arm sample celltype)
mkdir -p "$FT" "$B/logs"
cd "$B"
STAGE="${1:-submit}"

if [ "$STAGE" = "submit" ]; then
  q="qsub -terse -q slice"
  je=$($q -N fi_eval -l h_vmem=96G "$SELF" eval)
  jt=$($q -N fi_truth -l h_vmem=96G -hold_jid "$je" "$SELF" truth)
  jf=$($q -N fi_fisher -l h_vmem=32G -hold_jid "$jt" "$SELF" fisher)
  js=$($q -N fi_ucov_setup -hold_jid "$jf" "$SELF" ucov_setup)
  nb=$(ls $PM/02_gt_bams_umi_v6/merged_NTC_ctrl.bam $PM/02_gt_bams_umi_v6/T[1-4]_NTC_ctrl.bam \
          $PM/02_gt_bams_celltype_umi_v6/NTC_ctrl__*.bam | wc -l)
  ju=$($q -N fi_ucov -l h_vmem=4G -t 1-"$nb" -hold_jid "$js" "$SELF" ucov | cut -d. -f1)
  jc=$($q -N fi_ucov_collect -hold_jid "$ju" "$SELF" ucov_collect)
  jsc=$($q -N fi_score -l h_vmem=48G -l h_rt=12:00:00 -t 1-3 -hold_jid "$jc" "$SELF" score | cut -d. -f1)
  jr=$($q -N fi_rank -l h_vmem=16G -hold_jid "$jsc" "$SELF" rank)
  echo "eval $je | truth $jt | fisher $jf | ucov_setup $js | ucov $ju (1-$nb) | collect $jc | score $jsc | rank $jr"
  exit 0
fi

if [ "$STAGE" = "eval" ]; then GT_COVER_GATE=$GATE $RS human/scripts/gt_spectrum12_eval_v6.R; exit 0; fi

if [ "$STAGE" = "truth" ]; then
  [ -s "$GT_SP12_DIR/eval_$GATE/truth12_sites.tsv.gz" ] || { echo "ERROR: eval_$GATE missing" >&2; exit 1; }
  $RS human/scripts/gt_truth_agtc_test_v6.R
  $RS human/scripts/gt_truth_v7.R
  $RS human/scripts/gt_truth_v7_variants.R
  exit 0
fi

if [ "$STAGE" = "fisher" ]; then $RS human/scripts/gt_v8_fisher_truth.R; exit 0; fi

UC_ENV="OUT=$FT/unit_cov TRUTH=$TRUTH"
if [ "$STAGE" = "ucov_setup" ]; then env $UC_ENV bash "$PM/scripts/gt_v8_unit_cov.sh" setup; exit 0; fi
if [ "$STAGE" = "ucov" ]; then env $UC_ENV bash "$PM/scripts/gt_v8_unit_cov.sh"; exit 0; fi
if [ "$STAGE" = "ucov_collect" ]; then env $UC_ENV bash "$PM/scripts/gt_v8_unit_cov.sh" collect; exit 0; fi

if [ "$STAGE" = "score" ]; then
  tier=${TIERS[$(( SGE_TASK_ID - 1 ))]}
  GT_V7_TRUTH=$TRUTH GT_UNIT_COV=$FT/unit_cov/ko_truth_unit_cov.tsv.gz GT_UNIT_KEEP=${GT_UNIT_KEEP:-$V8/unit_cov/celltype_keep_min20.txt} \
    $RS human/scripts/gt_score_v7.R "$tier" "$PM/${PUB[$tier]}" "$FT/${SCN:-scoring}/mismatch/$tier" '_step5_nomm\.tsv\.gz$'
  exit 0
fi

if [ "$STAGE" = "rank" ]; then
  $RS human/scripts/gt_rank_v7.R "$FT/${SCN:-scoring}/mismatch"
  $RS human/scripts/gt_v8_fisher_compare.R
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
