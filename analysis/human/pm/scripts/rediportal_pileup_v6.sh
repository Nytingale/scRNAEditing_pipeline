#!/bin/bash
# rediportal_pileup_v6.sh   (HANDOFF 2026-09-29 UPDATE-25)
# Stranded coverage at every REDIportal v2 A-to-I position (chr1-22, X, Y) in each NTC unit
# BAM of the three GT tiers, for the REDIportal-as-truth benchmark (gt_score_rediportal_v6.R).
#   arm      02_gt_bams_umi_v6/merged_NTC_ctrl.bam
#   sample   02_gt_bams_umi_v6/T{1..4}_NTC_ctrl.bam
#   celltype 02_gt_bams_celltype_umi_v6/NTC_ctrl__*.bam (23)
# Same pileup as the truth build (samtools-1.20 mpileup -B, misc_scripts/parse_pileup_stranded.pl,
# libtype forward), BQ >= 30 and MQ >= 30 as the truth scan. The parser reports coverage on the transcribed strand
# implied by the reference base (A: plus-strand reads, T: minus-strand reads), i.e. REDIportal's strand.
#
#   bash rediportal_pileup_v6.sh setup                  # per-chromosome BEDs + task list
#   qsub -t 1-$(wc -l < $OUT/tasks.tsv) -tc 100 rediportal_pileup_v6.sh
# Output: gt_eval_v6/rediportal/pileup/<unit>/<chr>.tsv.gz  (chr pos ref coverage variant_cov)
# FULL=1 keeps the per-base columns (chr pos ref coverage variant_cov A C G T a c g t) in
# pileup_full/<unit>/ for the multiallelic count; merged NTC only = task lines 1-24:
#   qsub -t 1-24 -v FULL=1 rediportal_pileup_v6.sh
#$ -S /bin/bash
#$ -N redi_pileup
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=4G
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
OUT=$PM/gt_eval_v6/rediportal
VCF=$B/human/ref_annots/REDIportal_hg38_v2.vcf.gz
FASTA=$B/human/ref_annots/hg38.analysisSet.fa
PERL=$B/misc_scripts/parse_pileup_stranded.pl
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

if [ "${1:-}" = "setup" ]; then
  mkdir -p "$OUT/bed" "$OUT/pileup"
  zcat "$VCF" | awk -v OFS='\t' -v d="$OUT/bed" '!/^#/{print $1, $2-1, $2 > (d"/"$1".bed")}'
  for c in $(ls "$OUT/bed"); do case " $CHRS " in *" ${c%.bed} "*) ;; *) rm "$OUT/bed/$c";; esac; done
  : > "$OUT/tasks.tsv"
  for bam in "$PM/02_gt_bams_umi_v6/merged_NTC_ctrl.bam" "$PM"/02_gt_bams_umi_v6/T[1-4]_NTC_ctrl.bam \
             "$PM"/02_gt_bams_celltype_umi_v6/NTC_ctrl__*.bam; do
    for c in $CHRS; do printf '%s\t%s\n' "$bam" "$c" >> "$OUT/tasks.tsv"; done
  done
  echo "BED sites: $(cat "$OUT"/bed/*.bed | wc -l); tasks: $(wc -l < "$OUT/tasks.tsv")"
  exit 0
fi

read -r BAM CHR < <(sed -n "${SGE_TASK_ID}p" "$OUT/tasks.tsv")
unit=$(basename "$BAM" .bam)
if [ "${FULL:-0}" = "1" ]; then PD=pileup_full; COLS=1-13; else PD=pileup; COLS=1-5; fi
mkdir -p "$OUT/$PD/$unit"
o="$OUT/$PD/$unit/$CHR.tsv.gz"
samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$CHR" -l "$OUT/bed/$CHR.bed" -f "$FASTA" "$BAM" 2>/dev/null \
  | perl "$PERL" forward | cut -f"$COLS" | gzip > "$o.tmp"
mv "$o.tmp" "$o"
echo "$unit $CHR $(zcat "$o" | tail -n +2 | wc -l) covered positions"
