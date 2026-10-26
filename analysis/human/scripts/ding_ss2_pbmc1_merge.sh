#!/bin/bash
# ding_ss2_pbmc1_merge.sh   (2026-10-18; Reviewer 3, major 2)
# Merge the 384 PBMC1 Smart-seq2 cell BAMs (same inputs as human/jzlss2/wasp_refbias/merged_all.bam, PBMC1 only, so
# the donor matches the 10x 3' v3 PBMC1 data) and count primary alignments in it and in the 10x merged BAM.
#   qsub human/scripts/ding_ss2_pbmc1_merge.sh
#$ -S /bin/bash
#$ -N ding_p1m
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 4
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
set -euo pipefail
H=${SCRNA_DIR:-.}/human
W=$H/jzlss2/wasp_refbias
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
awk -F'\t' '$3 == "PBMC1" {print $1}' $H/jzlss2/name_match.txt | grep -F -f - $W/bamlist_all.txt > $W/bamlist_pbmc1.txt
test "$(wc -l < $W/bamlist_pbmc1.txt)" -eq 384
samtools merge -f -@ 4 -b $W/bamlist_pbmc1.txt $W/merged_pbmc1.bam
samtools index -@ 4 $W/merged_pbmc1.bam
{
  echo -e "bam\treads_primary"
  for b in $H/jzlxv3/wasp_refbias/merged_all.bam $W/merged_all.bam $W/merged_pbmc1.bam; do
    echo -e "$b\t$(samtools view -c -@ 4 -F 0x904 $b)"
  done
} > $H/jzlxv3/three_prime_bias/read_counts.tsv
cat $H/jzlxv3/three_prime_bias/read_counts.tsv
