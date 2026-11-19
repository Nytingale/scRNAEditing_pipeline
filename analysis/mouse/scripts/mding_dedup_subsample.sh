#!/bin/bash
# mding_dedup_subsample.sh   (2026-11-07; Reviewer 3, major 2; mouse version of human/scripts/ding_ss2_dedup_subsample.sh)
# Ding mouse cortex, TSO-clipped runs: equal unique reads on both platforms per mouse, from the pipeline's deduplicated
# BAMs: Smart-seq2 mouse/mjzlss2/results/<M>/02_dedup/<M>_rmDup.bam (per-cell coordinate deduplication; paired-end) and
# 10x 3' mouse/mjzl10x/results/<M>/02_dedup/<M>_umiMd.bam (UMI-aware MarkDuplicates; single-end; 2026-11-14, was the
# umi_tools <M>_umiDup.bam with outputs in 02_dedup/subsample, kept).
# K = the smaller of (10x primary mapped reads, 2 x Smart-seq2 proper pairs), made even; both brought to exactly K reads:
#   Smart-seq2 K / 2 pairs sampled as whole pairs, seeds 20261024/25/26; 10x K reads, seed 20261024
# Selection sampling (Knuth, Algorithm S); output stays coordinate-sorted.
#   tasks 1-8: mouse (Cortex1: 1-4, Cortex2: 5-8) x (Smart-seq2 seed 1-3, 10x)
#   qsub -t 1-8 mouse/scripts/mding_dedup_subsample.sh
# Output: mouse/mjzlss2/results/<M>/02_dedup/subsample_umimd/<M>_rmDup_sub<seed>.bam
#         mouse/mjzl10x/results/<M>/02_dedup/subsample_umimd/<M>_umiMd_sub<seed>.bam   (+ .bai, .count)
#$ -S /bin/bash
#$ -N mding_sub
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/mjzl10x/three_prime_bias/logs/
#$ -j y
#$ -pe threaded 4
#$ -l h_vmem=4G
#$ -l h_rt=24:00:00
set -euo pipefail
MO=${SCRNA_DIR:-.}/mouse
MICE=(Cortex1 Cortex2)
SEEDS=(20261024 20261025 20261026 20261024)
i=$((SGE_TASK_ID - 1)); M=${MICE[$((i / 4))]}; j=$((i % 4)); SEED=${SEEDS[$j]}
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
SS=$MO/mjzlss2/results/$M/02_dedup/${M}_rmDup.bam
TX=$MO/mjzl10x/results/$M/02_dedup/${M}_umiMd.bam
N10=$(samtools view -c -@ 4 -F 0x904 "$TX")
NSP=$(samtools view -c -@ 4 -f 0x43 -F 0x90C "$SS")
K=$(( (N10 < 2 * NSP ? N10 : 2 * NSP) / 2 * 2 ))
echo "=== $(date) $M: 10x $N10 reads, Smart-seq2 $NSP proper pairs -> K = $K"

if [ "$j" -le 2 ]; then
  IN=$SS; OUT=$MO/mjzlss2/results/$M/02_dedup/subsample_umimd; o=$OUT/${M}_rmDup_sub${SEED}.bam
  FILT="-f 0x3 -F 0x90C"; N=$NSP
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
  IN=$TX; OUT=$MO/mjzl10x/results/$M/02_dedup/subsample_umimd; o=$OUT/${M}_umiMd_sub${SEED}.bam
  FILT="-F 0x904"; N=$N10
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
echo -e "seed\treads_primary\tsource_units\tK\n$SEED\t$n\t$N\t$K" > "${o%.bam}.count"
[ "$n" -eq "$K" ] || { echo "count $n != $K"; exit 1; }
echo "=== $(date) done $n reads"
