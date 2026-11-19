#!/bin/bash
# ding_ss2_dedup_subsample.sh   (2026-10-24; Reviewer 3, major 2; 10x UMI-aware MarkDuplicates 2026-11-14)
# Equal unique reads on both platforms without per-site thinning, from the pipeline's own deduplicated BAMs:
# Smart-seq2 PBMC1 results_v8/PBMC1/02_dedup/PBMC1_rmDup.bam (GATK MarkDuplicates per cell, duplicates removed;
# paired-end) and 10x 3' v3 02_bams_umimd/PBMC1.bam (UMI-aware MarkDuplicates per cell, pbmc_umimd.sh; single-end,
# 217,323,582 primary mapped reads). Both are brought to exactly K reads (the 10x count made even; qsub -v K):
#   tasks 1-3  Smart-seq2: proper pairs only (samtools -f 0x3 -F 0x90C), K / 2 pairs sampled as whole pairs (both
#              mates kept), seeds 20261024/25/26
#   task 4     10x: K of its primary mapped reads (-F 0x904; all of them when the count is even), seed 20261024
# (Until 2026-11-14 the 10x input was the umi_tools PBMC1_umiDup.bam and the outputs went to 02_dedup/subsample.)
# Selection sampling (Knuth, Algorithm S) over pairs (decided at the first mate, applied to the second) or reads gives an
# exact, uniformly random subset; output stays coordinate-sorted.
#   qsub -t 1-4 human/scripts/ding_ss2_dedup_subsample.sh
# Output: human/jzlss2/results_v8/PBMC1/02_dedup/subsample_umimd/PBMC1_rmDup_sub<seed>.bam
#         human/jzlxv3/02_bams_umimd/subsample/PBMC1_umiMd_sub<seed>.bam   (+ .bai, .count)
#$ -S /bin/bash
#$ -N ding_sub
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 4
#$ -l h_vmem=4G
#$ -l h_rt=24:00:00
set -euo pipefail
H=${SCRNA_DIR:-.}/human
# K = the 10x UMI-aware MarkDuplicates primary mapped count made even (217323582); qsub -v K=<K>
K=${K:?K = even 10x UMI-aware MarkDuplicates primary mapped read count}
SEEDS=(20261024 20261025 20261026 20261024)
SEED=${SEEDS[$((SGE_TASK_ID - 1))]}
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"

if [ "$SGE_TASK_ID" -le 3 ]; then
  IN=$H/jzlss2/results_v8/PBMC1/02_dedup/PBMC1_rmDup.bam
  OUT=$H/jzlss2/results_v8/PBMC1/02_dedup/subsample_umimd; o=$OUT/PBMC1_rmDup_sub${SEED}.bam
  FILT="-f 0x3 -F 0x90C"
  N=$(samtools view -c -@ 4 -f 0x43 -F 0x90C "$IN")
  echo "=== $(date) Smart-seq2 seed $SEED: $((K / 2)) of $N proper pairs"
  mkdir -p "$OUT"
  samtools view -h -@ 2 $FILT "$IN" \
    | awk -v N=$N -v K=$((K / 2)) -v S=$SEED 'BEGIN { srand(S) }
        /^@/ { print; next }
        $1 in keep { if (keep[$1]) print; delete keep[$1]; next }
        { k = (rand() * (N - seen) < K - sel); keep[$1] = k; sel += k; seen++; if (k) print }
        END { n = 0; for (q in keep) n++
              if (seen != N || sel != K || n) { print "pairs " seen " selected " sel " unmatched " n > "/dev/stderr"; exit 1 } }' \
    | samtools view -b -@ 2 -o "$o.tmp" -
else
  IN=$H/jzlxv3/02_bams_umimd/PBMC1.bam
  OUT=$H/jzlxv3/02_bams_umimd/subsample; o=$OUT/PBMC1_umiMd_sub${SEED}.bam
  FILT="-F 0x904"
  N=$(samtools view -c -@ 4 $FILT "$IN")
  echo "=== $(date) 10x seed $SEED: $K of $N reads"
  mkdir -p "$OUT"
  samtools view -h -@ 2 $FILT "$IN" \
    | awk -v N=$N -v K=$K -v S=$SEED 'BEGIN { srand(S) }
        /^@/ { print; next }
        { if (rand() * (N - seen) < K - sel) { print; sel++ }; seen++ }
        END { if (seen != N || sel != K) { print "seen " seen " selected " sel > "/dev/stderr"; exit 1 } }' \
    | samtools view -b -@ 2 -o "$o.tmp" -
fi
mv "$o.tmp" "$o"
samtools index -@ 4 "$o"
n=$(samtools view -c -@ 4 $FILT "$o")
echo -e "seed\treads_primary\tsource_units\n$SEED\t$n\t$N" > "${o%.bam}.count"
[ "$n" -eq "$K" ] || { echo "count $n != $K"; exit 1; }
echo "=== $(date) done $n reads"
