#!/bin/bash
# ding_dedup_read_counts.sh   (2026-10-19; Reviewer 3, major 2)
# Primary mapped reads (samtools view -c -F 0x904, as ding_ss2_pbmc1_merge.sh) of the deduplicated Ding PBMC1 BAMs
# (10x UMI-deduplicated, Smart-seq2 per-cell coordinate-deduplicated), appended to three_prime_bias/read_counts.tsv
# for the thinning in ding_3prime_bias_denovo.R (X10_PLATFORM=x10u SS2_PLATFORM=ss2p1d).
#   qsub human/scripts/ding_dedup_read_counts.sh
#$ -S /bin/bash
#$ -N ding_rc
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 4
#$ -l h_vmem=2G
#$ -l h_rt=4:00:00
set -euo pipefail
H=${SCRNA_DIR:-.}/human
RC=$H/jzlxv3/three_prime_bias/read_counts.tsv
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
for b in $H/jzlxv3/results_v8/PBMC1/02_dedup/PBMC1_umiDup.bam $H/jzlss2/results_v8/PBMC1/02_dedup/PBMC1_rmDup.bam; do
  grep -qF "$b" "$RC" && continue
  n=$(samtools view -c -@ 4 -F 0x904 "$b")
  echo -e "$b\t$n" >> "$RC"
done
cat "$RC"
