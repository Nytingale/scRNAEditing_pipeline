#!/bin/bash
# =============================================================================
# gt_v9_m_null.sh  [submit|eval|truth]   (2026-11-16; Response Figure 1b on the v9 truth and the seed-8787 pileups)
# Site universe of the e^-m null check (GTDEF_v8_m_null_check_plot.R with GT_V9=1): the NTC rules with the gate
# relaxed to NTC cov >= 10 and KO cov >= 1 (gate depth1), as run_gt_v8_fisher.sh eval + truth, with the v9 site
# definition (GT_V7_VARIANT=AGTC: read-orientation A>G or T>C, one orientation per position, then multiallelic,
# PGP1 and nearby-mismatch cleaning with w1 searched).
#   eval   gt_spectrum12_eval_v6.R GT_COVER_GATE=depth1 -> gt_eval_v6/spectrum12/eval_depth1/
#   truth  gt_truth_agtc_test_v6.R -> gt_truth_v7.R -> gt_truth_v7_variants.R -> gt_eval_v9/m_null/
# Usage : bash human/pm/scripts/gt_v9_m_null.sh submit   (login node; qsub shim on PATH)
#$ -S /bin/bash
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -l h_vmem=96G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
MN=$PM/gt_eval_v9/m_null
SELF=$PM/scripts/gt_v9_m_null.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH="${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:${PATH}"
GATE=depth1
export GT_SP12_DIR=$PM/gt_eval_v6/spectrum12 GT_AGTC_DIR=$MN/agtc_test GT_V7_TRUTH_DIR=$MN/truth_v7 \
       GT_V7_VARIANTS_DIR=$MN/truth GT_GATES=$GATE GT_TRUTH_ONLY=1 GT_V7_VARIANT=AGTC
mkdir -p "$MN" "$B/logs"
cd "$B"
STAGE="${1:-submit}"

if [ "$STAGE" = "submit" ]; then
  je=$(qsub -terse -N mn9_eval "$SELF" eval)
  jt=$(qsub -terse -N mn9_truth -hold_jid "$je" "$SELF" truth)
  echo "eval $je | truth $jt"
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
