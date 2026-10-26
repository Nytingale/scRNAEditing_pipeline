#!/bin/bash
# gt_v8_unit_cov.sh   (2026-10-13)
# Stranded coverage at every v8 KO truth position (gt_eval_v8/truth/gt_v7_truth_sites.tsv.gz) in each NTC
# unit BAM of the three GT tiers, for the per-unit recall denominator (gt_score_v7.R, GT_UNIT_COV).
# Same pileup as the truth build and rediportal_pileup_v6.sh (samtools-1.20 mpileup -B, BQ >= 30, MQ >= 30,
# misc_scripts/parse_pileup_stranded.pl forward: coverage on the strand implied by the reference base).
#   arm      02_gt_bams_umi_v6/merged_NTC_ctrl.bam      -> unit NTC_ctrl
#   sample   02_gt_bams_umi_v6/T{1..4}_NTC_ctrl.bam     -> unit T<n>_NTC_ctrl
#   celltype 02_gt_bams_celltype_umi_v6/NTC_ctrl__*.bam -> unit NTC_ctrl__<celltype>
# (REDIportal positions: gt_eval_v6/rediportal/pileup/<unit>/, already per unit.)
#
#   bash gt_v8_unit_cov.sh setup
#   qsub -t 1-28 gt_v8_unit_cov.sh
#   bash gt_v8_unit_cov.sh collect      # -> $OUT/ko_truth_unit_cov.tsv.gz (unit site coverage variant_cov)
#$ -S /bin/bash
#$ -N gt_v8_ucov
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=4G
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
# split-by-cells test: OUT=split_cells/unit_cov TRUTH=split_cells/truth/v7_nomm_pgp1/... BAMS="split_cells/bams/umi/*_B.bam"
# (unit names drop the _B suffix)
OUT=${OUT:-$PM/gt_eval_v8/unit_cov}
TRUTH=${TRUTH:-$PM/gt_eval_v8/truth/gt_v7_truth_sites.tsv.gz}
BAMS=${BAMS:-"$PM/02_gt_bams_umi_v6/merged_NTC_ctrl.bam $PM/02_gt_bams_umi_v6/T[1-4]_NTC_ctrl.bam $PM/02_gt_bams_celltype_umi_v6/NTC_ctrl__*.bam"}
FASTA=$B/human/ref_annots/hg38.analysisSet.fa
PERL=$B/misc_scripts/parse_pileup_stranded.pl
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"

if [ "${1:-}" = "setup" ]; then
  mkdir -p "$OUT/ko_truth"
  for t in $TRUTH; do zcat "$t" | awk -v OFS='\t' 'NR > 1 {print $3, $4 - 1, $4}'; done \
    | sort -k1,1 -k2,2n -u > "$OUT/ko_truth_sites.bed"
  ls $BAMS > "$OUT/tasks.tsv"
  echo "positions: $(wc -l < "$OUT/ko_truth_sites.bed"); BAMs: $(wc -l < "$OUT/tasks.tsv")"
  exit 0
fi

if [ "${1:-}" = "collect" ]; then
  { printf 'unit\tsite\tcoverage\tvariant_cov\n'
    for f in "$OUT"/ko_truth/*.tsv.gz; do
      u=$(basename "$f" .tsv.gz); u=${u%_B}; [ "$u" = merged_NTC_ctrl ] && u=NTC_ctrl
      zcat "$f" | awk -v OFS='\t' -v u="$u" 'NR > 1 {print u, $1":"$2, $4, $5}'
    done; } | gzip > "$OUT/ko_truth_unit_cov.tsv.gz"
  zcat "$OUT/ko_truth_unit_cov.tsv.gz" | awk 'NR > 1 {n[$1]++; if ($3 >= 10) c[$1]++} END {for (u in n) print u, n[u], c[u] + 0}' | sort
  exit 0
fi

BAM=$(sed -n "${SGE_TASK_ID}p" "$OUT/tasks.tsv")
unit=$(basename "$BAM" .bam)
o="$OUT/ko_truth/$unit.tsv.gz"
samtools mpileup -B -q 30 -Q 30 -d 1000000 -l "$OUT/ko_truth_sites.bed" -f "$FASTA" "$BAM" 2>/dev/null \
  | perl "$PERL" forward | cut -f1-5 | gzip > "$o.tmp"
mv "$o.tmp" "$o"
echo "$unit $(zcat "$o" | tail -n +2 | wc -l) covered positions"
