#!/bin/bash
# ding_3prime_bias_pileup_sub.sh   (2026-10-24; Reviewer 3, major 2)
# ding_3prime_bias_pileup.sh on the subsampled deduplicated BAMs (ding_ss2_dedup_subsample.sh), same mpileup and
# parsing: 3 Smart-seq2 PBMC1 seeds (platforms ss2p1dsm<seed>) and the 10x UMI-aware MarkDuplicates one (x10ms<seed>).
# (Until 2026-11-14: ss2p1ds<seed> and the umi_tools x10us<seed>; those pileups are kept.)
#   qsub -t 1-92 -hold_jid ding_sub human/scripts/ding_3prime_bias_pileup_sub.sh   (23 chromosomes x 4 BAMs)
#$ -S /bin/bash
#$ -N ding_pus
# -q <queue>: add as a "#$" line if your cluster needs one
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

CHRS=($(for i in $(seq 1 22) X; do echo chr$i; done))
SEEDS=(20261024 20261025 20261026)
i=$((SGE_TASK_ID - 1))
j=$((i / 23)); c=${CHRS[$((i % 23))]}
if [ "$j" -le 2 ]; then
  s=${SEEDS[$j]}; bam=$H/jzlss2/results_v8/PBMC1/02_dedup/subsample_umimd/PBMC1_rmDup_sub${s}.bam; plat=ss2p1dsm$s
else
  s=${SEEDS[0]}; bam=$H/jzlxv3/02_bams_umimd/subsample/PBMC1_umiMd_sub${s}.bam; plat=x10ms$s
fi
o="$OUT/${plat}_${c}.tsv.gz"
echo "=== $(date) $bam $c"
samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$c" -f "$FASTA" "$bam" 2>/dev/null \
  | awk -F'\t' '$4 >= 10 && ($3 == "A" || $3 == "T" || $3 == "a" || $3 == "t")' \
  | perl "$PERL" unstranded | awk -F'\t' -v OFS='\t' 'NR == 1 || $4 >= 10 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' \
  | gzip > "$o.tmp"
mv "$o.tmp" "$o"
echo "=== $(date) done $(zcat "$o" | tail -n +2 | wc -l) positions"
