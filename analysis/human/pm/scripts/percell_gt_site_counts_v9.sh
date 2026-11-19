#!/bin/bash
# v9 copy of percell_gt_site_counts_v8.sh (2026-11-17): sites and classes from gt_eval_v9/truth/gt_v9_positive_sites.tsv.gz
# (STAGE=sites below; union over gates x t25/t50/t75). Per-cell ref/edited read counts at the KO ground-truth
# positives, one array task per per-teratoma GT sample (T1-T4 x 5 arms = 20). Unprocessed reads (02_gt_bams), as
# in the v5 run (percell_gt_site_counts.sh). Per-truth membership: percell_gt_coverage_v9/gt_v9_site_classes.tsv.
# Usage:  STAGE=sites bash human/pm/scripts/percell_gt_site_counts_v9.sh; qsub -t 1-20 human/pm/scripts/percell_gt_site_counts_v9.sh
#$ -N pc_gtcov9
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 1
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
# Do NOT pass -V (SGE JSV / truncated PATH).

set -eo pipefail
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH

BASE="${SCRNA_DIR:-.}"
PM="${BASE}/human/pm"
OUT="${PM}/percell_gt_coverage_v9"
SITES="${OUT}/gt_sites_v9.tsv"
if [ "${STAGE:-}" = sites ]; then
  mkdir -p "$OUT"
  zcat "$PM/gt_eval_v9/truth/gt_v9_positive_sites.tsv.gz" |
    awk -F'\t' -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "gate","site","class_t25","class_t50","class_t75"; next}
      {print $h["gate"],$h["site"],$h["class_t25"],$h["class_t50"],$h["class_t75"]}' > "$OUT/gt_v9_site_classes.tsv"
  { printf "site\tclass\n"; awk -F'\t' 'NR>1 && ($3=="positive"||$4=="positive"||$5=="positive"){print $2"\tTRUE"}' \
      "$OUT/gt_v9_site_classes.tsv" | sort -u; } > "$SITES"
  echo "$(( $(wc -l < "$SITES") - 1 )) union positives -> $SITES"; exit 0
fi
[ -s "$SITES" ] || { echo "missing $SITES" >&2; exit 1; }

SAMPLES=($(for t in T1 T2 T3 T4; do for a in AAVS1_ctrl ADAR1_KO ADARB1_KO ADARB2_KO NTC_ctrl; do echo "${t}_${a}"; done; done))
S="${SAMPLES[$((SGE_TASK_ID - 1))]}"
BAM="${PM}/02_gt_bams/${S}.bam"
echo "=== $(date) === ${S} (${BAM})"

python "${PM}/scripts/percell_gt_site_counts.py" \
    --bam "$BAM" \
    --sample "$S" \
    --cells "${PM}/gt_barcodes/${S}.tsv" \
    --sites "$SITES" \
    --fasta "${BASE}/human/ref_annots/hg38.analysisSet.fa" \
    --out "${OUT}/${S}_gt_site_counts.tsv.gz"

echo "=== $(date) === done ${S}"
