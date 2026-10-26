#!/bin/bash
# run_gt_rediportal_v7.sh   (2026-10-08)
# Score the published teratoma call sets (RNA-frame-fixed REDItools3) against the three
# REDIportal positive variants used in the nine-truth suite, writing under gt_eval_v7.
# Positives and pileups are reused from gt_eval_v6/rediportal/ (same merged-NTC coverage
# gates as gt_score_rediportal_v6.R: cov10_alt1 / alt3 / alt3_e05). The v7 clean figures read
# scoring_mismatch (run_gt_rediportal_v7_states.sh mismatch), not this output.
#   qsub human/pm/scripts/run_gt_rediportal_v7.sh
#$ -S /bin/bash
#$ -N gt_redi_v7
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
V6=$B/human/pm/gt_eval_v6/rediportal
OUT=$B/human/pm/gt_eval_v7/rediportal/scoring
cd "$B"
# pileups + dbSNP-common flag from the v6 REDIportal build (after_rediportal_pileup_v6.sh)
[ -s "$V6/rediportal_dbsnp_common_sites.txt" ] || {
  echo "missing $V6/rediportal_dbsnp_common_sites.txt; run after_rediportal_pileup_v6.sh first" >&2
  exit 1
}
n_pu=$(find "$V6/pileup" -name 'chr*.tsv.gz' | wc -l)
[ "$n_pu" -eq "$(wc -l < "$V6/tasks.tsv")" ] || {
  echo "ERROR: $n_pu pileup outputs, expected $(wc -l < "$V6/tasks.tsv")" >&2
  exit 1
}
mkdir -p "$OUT"
export GT_REDI_OUT="$OUT"
# PGP1 SNV positions removed by default (as run_gt_rediportal_v7_states.sh); multiallelic sites
# removed by gt_score_rediportal_v6.R default
export GT_REDI_EXCLUDE_BED=${GT_REDI_EXCLUDE_BED-$B/human/pm/pgp1_snps/PGP1_hu43860C_hg38_snvs.bed}
echo "GT_REDI_OUT=$GT_REDI_OUT GT_REDI_EXCLUDE_BED=$GT_REDI_EXCLUDE_BED GT_REDI_DROP_MULTI=${GT_REDI_DROP_MULTI:-1} (call trees = published results_*_v6)"
$RS human/scripts/gt_score_rediportal_v6.R
for v in cov10_alt1 alt3 alt3_e05; do
  f=$OUT/gt_tier_scoring_$v.tsv
  [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
  echo "$v: $(wc -l < "$f") rows"
done
echo "done"
