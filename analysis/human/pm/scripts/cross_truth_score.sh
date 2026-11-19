#!/bin/bash
# cross_truth_score.sh
# Scores every call set against every truth, to measure how much the deduplication used to build the truth favours
# the configurations called on the same BAMs (shared-read bias).
#   truths     v9 (published, unseeded umi_tools BAMs), seed1, seed2 (umi_tools --random-seed 1 / 2), umimd
#              (UMI-aware MarkDuplicates BAMs); built by umi_seed_reps.sh
#   call sets  pub (58 published configurations), ch (2 Cellsnp-lite hybrid), umimd (16), seed1, seed2 (4 BCFtools
#              umi_tools configurations re-called on the seed BAMs)
#   tiers      teratoma arm, sample, celltype; aorta arm, sample
# Same settings as gt_v9_build.sh score (teratoma) and aorta_v9_build.sh score (aorta).
# Array task -> (dataset, truth, call set, tier): tasks 1-60 teratoma, 61-100 aorta; 101-125 agreement truth
#   (cross_truth_agreement.R; STAGE=submit_agree builds it, then scores 101-115 teratoma, 116-125 aorta);
#   126-150 coordinate truth (umi_seed_reps.sh SEED=coord; STAGE=submit_coord HOLD=<truth jobs>).
# Output: human/pm/gt_eval_v9/cross_truth/scores/<ds>/<truth>/<calls>/<tier>/gt_v7_scores.tsv; then
#         Rscript human/scripts/cross_truth_summary.R
# Usage:  STAGE=submit bash human/pm/scripts/cross_truth_score.sh   (login node: qsub only; HOLD=<jobids> optional)
#$ -S /bin/bash
#$ -N cross_truth
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/human/pm/gt_eval_v9/cross_truth/logs/
#$ -j y
#$ -l h_vmem=16G
#$ -l h_rt=4:00:00
set -euo pipefail
export LC_ALL=C
B=${SCRNA_DIR:-.}
PM=$B/human/pm
A=$B/mouse/aorta_gse255661
V8=$PM/gt_eval_v8
O=$B/human/pm/gt_eval_v9/cross_truth
SELF=$PM/scripts/cross_truth_score.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
TRUTHS=(v9 seed1 seed2 umimd)
CALLS=(pub ch umimd seed1 seed2)
STAGE=${STAGE:-${1:-score}}
cd "$B"

Q=$(command -v qsub || echo /cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub)
case $STAGE in
submit)
  mkdir -p $O/logs
  $Q -terse -t 1-100 -tc 40 ${HOLD:+-hold_jid $HOLD} $SELF
  exit 0 ;;
submit_coord)
  mkdir -p $O/logs
  $Q -terse -t 126-150 ${HOLD:+-hold_jid $HOLD} $SELF
  exit 0 ;;
submit_agree)
  mkdir -p $O/logs
  j=$($Q -terse -v STAGE=agree -t 1-2 -l h_vmem=32G $SELF | cut -d. -f1)
  $Q -terse -t 101-125 -hold_jid $j $SELF
  exit 0 ;;
agree)
  ds=$( [ "$SGE_TASK_ID" = 1 ] && echo teratoma || echo aorta )
  $RS human/scripts/cross_truth_agreement.R $ds
  exit 0 ;;
esac

i=$(( SGE_TASK_ID - 1 ))
if [ $i -ge 100 ]; then
  # tasks 101-125: agreement truth (cross_truth_agreement.R); 126-150: coordinate truth (umi_seed_reps.sh SEED=coord);
  # within each block 15 teratoma, then 10 aorta
  TRUTHS=(agree); [ $i -ge 125 ] && TRUTHS=(coord)
  i=$(( (i - 100) % 25 ))
  if [ $i -lt 15 ]; then ds=teratoma; TIERS=(arm sample celltype); else ds=aorta; TIERS=(arm sample); i=$(( i - 15 )); fi
elif [ $i -lt 60 ]; then
  ds=teratoma; TIERS=(arm sample celltype)
else
  ds=aorta; TIERS=(arm sample); i=$(( i - 60 ))
fi
nt=${#TIERS[@]}
tier=${TIERS[$(( i % nt ))]}
calls=${CALLS[$(( i / nt % 5 ))]}
truth=${TRUTHS[$(( i / nt / 5 ))]}
out=$O/scores/$ds/$truth/$calls/$tier
mkdir -p "$out"
echo "task $SGE_TASK_ID: $ds truth=$truth calls=$calls tier=$tier -> $out"

if [ $ds = teratoma ]; then
  declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
  case $calls in
    pub)   res=$PM/${PUB[$tier]} ;;
    ch)    res=$V8/cellsnp_hybrid/results/centre/${PUB[$tier]} ;;
    umimd) res=$PM/results_umimd_$tier ;;
    seed*) res=$PM/umi_seed_reps/$calls/results_$tier ;;
  esac
  T=$( [ $truth = v9 ] && echo $PM/gt_eval_v9 || echo $PM/umi_seed_reps/$truth/v9 )
  export GT_UNIT_RULE=gate_cov GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pwcov.txt GT_FP_RULE=ag
  export GT_V7_TRUTH=$T/truth/gt_v9_truth_sites.tsv.gz GT_UNIT_COV=$T/unit_cov/ko_truth_unit_cov.tsv.gz
  if [ $truth = agree ]; then
    export GT_V7_TRUTH=$O/truth_agree/$ds/truth_sites.tsv.gz GT_UNIT_COV=$O/truth_agree/$ds/unit_cov.tsv.gz
  fi
else
  case $calls in
    pub)   res=$A/results_v8 ;;
    ch)    res=$A/gt_eval_v8/cellsnp_hybrid/results/results_v8 ;;
    umimd) res=$A/results_v8_umimd ;;
    seed*) res=$A/umi_seed_reps/$calls/results_aorta ;;
  esac
  T=$( [ $truth = v9 ] && echo $A/truth/eval_v9/excluded || echo $A/umi_seed_reps/$truth/eval_v9/excluded )
  export GT_FP_RULE=ag GT_UNIT_RULE=gate_cov GT_UNIT_EDIT_UNIT=CTRL
  export GT_UNIT_RE=$( [ "$tier" = arm ] && echo '^CTRL$' || echo '^Control_' )
  export GT_V7_TRUTH=$T/gt_aorta_truth_sites.tsv.gz GT_UNIT_COV=$T/unit_cov.tsv.gz
  if [ $truth = agree ]; then
    export GT_V7_TRUTH=$O/truth_agree/$ds/truth_sites.tsv.gz GT_UNIT_COV=$O/truth_agree/$ds/unit_cov.tsv.gz
  fi
fi
for f in "$GT_V7_TRUTH" "$GT_UNIT_COV"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done
n=$(find "$res" -name '*_step5_nomm.tsv.gz' | wc -l)
[ "$n" -gt 0 ] || { echo "no step5_nomm tables under $res" >&2; exit 1; }
echo "truth $GT_V7_TRUTH; $n step5_nomm tables under $res"
$RS human/scripts/gt_score_v7.R "$tier" "$res" "$out" '_step5_nomm\.tsv\.gz$'
