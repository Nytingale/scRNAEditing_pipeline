#!/bin/bash
# ding_3prime_bias_pileup.sh   (2026-10-18; Reviewer 3, major 2)
# Genome-wide pileup of Ding et al. human PBMC 10x 3' v3 (PBMC1) and Smart-seq2 (768 cells) merged BAMs (the WASP
# analysis inputs, all reads), for de novo A-to-G candidate sites beyond REDIportal (ding_3prime_bias_denovo.R).
# samtools-1.20 mpileup -B -q 30 -Q 30, parse_pileup_stranded.pl unstranded (Smart-seq2 is unstranded; both
# platforms parsed the same way), A/T reference positions with >= 10 reads kept.
#   qsub -t 1-115 human/scripts/ding_3prime_bias_pileup.sh     (23 chromosomes x 5 platforms; ss2p1 = Smart-seq2 PBMC1 only;
#   x10u = 10x UMI-deduplicated, ss2p1d = Smart-seq2 PBMC1 per-cell coordinate-deduplicated, from run_pipeline_dataset_v8.sh)
#$ -S /bin/bash
#$ -N ding_pu
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
H=$B/human
OUT=$H/jzlxv3/three_prime_bias/pileup
FASTA=$H/ref_annots/hg38.analysisSet.fa
PERL=$B/misc_scripts/parse_pileup_stranded.pl
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
mkdir -p "$OUT"

CHRS=($(for i in $(seq 1 22) X; do echo chr$i; done))
PLAT=(x10 ss2 ss2p1 x10u ss2p1d)
BAMS=($H/jzlxv3/wasp_refbias/merged_all.bam $H/jzlss2/wasp_refbias/merged_all.bam $H/jzlss2/wasp_refbias/merged_pbmc1.bam $H/jzlxv3/results_v8/PBMC1/02_dedup/PBMC1_umiDup.bam $H/jzlss2/results_v8/PBMC1/02_dedup/PBMC1_rmDup.bam)
i=$((SGE_TASK_ID - 1))
p=$((i / 23)); c=${CHRS[$((i % 23))]}
o="$OUT/${PLAT[$p]}_${c}.tsv.gz"
echo "=== $(date) ${PLAT[$p]} $c"
samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$c" -f "$FASTA" "${BAMS[$p]}" 2>/dev/null \
  | awk -F'\t' '$4 >= 10 && ($3 == "A" || $3 == "T" || $3 == "a" || $3 == "t")' \
  | perl "$PERL" unstranded | awk -F'\t' -v OFS='\t' 'NR == 1 || $4 >= 10 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' \
  | gzip > "$o.tmp"
mv "$o.tmp" "$o"
echo "=== $(date) done $(zcat "$o" | tail -n +2 | wc -l) positions"
