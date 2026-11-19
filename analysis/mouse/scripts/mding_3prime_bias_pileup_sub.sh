#!/bin/bash
# mding_3prime_bias_pileup_sub.sh   (2026-11-07; mouse version of human/scripts/ding_3prime_bias_pileup_sub.sh)
# Genome-wide unstranded pileup (A/T reference, >= 10 reads; samtools mpileup -B, BQ/MQ >= 30) of the subsampled
# deduplicated BAMs (mding_dedup_subsample.sh): per mouse, 3 Smart-seq2 seeds (ss2dsm<seed>) and the 10x UMI-aware
# MarkDuplicates one (x10ms<seed>) (until 2026-11-14: ss2ds<seed> / umi_tools x10us<seed>; those pileups are kept).
#   tasks 1-160: mouse (2) x BAM (4) x chromosome (chr1-19, chrX)
#   qsub -t 1-160 -hold_jid mding_sub mouse/scripts/mding_3prime_bias_pileup_sub.sh
# Output: mouse/mjzl10x/three_prime_bias/pileup/<M>_<plat>_<chr>.tsv.gz
#$ -S /bin/bash
#$ -N mding_pus
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/mjzl10x/three_prime_bias/logs/
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
MO=$B/mouse
OUT=$MO/mjzl10x/three_prime_bias/pileup
FASTA=$MO/ref_annots/mm39.fa
PERL=$B/misc_scripts/parse_pileup_stranded.pl
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
mkdir -p "$OUT"

MICE=(Cortex1 Cortex2)
CHRS=($(for i in $(seq 1 19) X; do echo chr$i; done))
SEEDS=(20261024 20261025 20261026)
i=$((SGE_TASK_ID - 1))
M=${MICE[$((i / 80))]}; r=$((i % 80)); j=$((r / 20)); c=${CHRS[$((r % 20))]}
if [ "$j" -le 2 ]; then
  s=${SEEDS[$j]}; bam=$MO/mjzlss2/results/$M/02_dedup/subsample_umimd/${M}_rmDup_sub${s}.bam; plat=ss2dsm$s
else
  s=${SEEDS[0]}; bam=$MO/mjzl10x/results/$M/02_dedup/subsample_umimd/${M}_umiMd_sub${s}.bam; plat=x10ms$s
fi
o="$OUT/${M}_${plat}_${c}.tsv.gz"
echo "=== $(date) $bam $c"
samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$c" -f "$FASTA" "$bam" 2>/dev/null \
  | awk -F'\t' '$4 >= 10 && ($3 == "A" || $3 == "T" || $3 == "a" || $3 == "t")' \
  | perl "$PERL" unstranded | awk -F'\t' -v OFS='\t' 'NR == 1 || $4 >= 10 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' \
  | gzip > "$o.tmp"
mv "$o.tmp" "$o"
echo "=== $(date) done $(zcat "$o" | tail -n +2 | wc -l) positions"
