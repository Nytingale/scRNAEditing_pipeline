#!/bin/bash
# =============================================================================
# merge_gt_umi_bams_v6.sh   (SGE array 1-5, one task per CRISPR arm)
#
# Rebuild the arm-level UMI BAMs as "dedup per teratoma (GEM well), THEN merge".
#
# Why: 02_gt_bams_umi/merged_<ARM>.bam (used by v5) is UMI_DEDUP run AFTER
# pooling T1-T4 (merge_gt_bams.sh -> results_gt_v2/<ARM>/02_dedup). Raw CR
# barcodes are not unique across teratomas, so reads from different teratomas
# sharing CR+UB+position could be collapsed (measured: 27 of 208M reads).
# Each teratoma is one GEM well (its two SRR runs share 43-72% of top-3000
# barcodes), so the per-teratoma UMI BAMs are the correct dedup unit.
#
# Inputs : human/pm/02_gt_bams_umi/T{1-4}_<ARM>.bam  (-> results_gt_persample/*/02_dedup/*_umiDup.bam)
# Outputs: human/pm/02_gt_bams_umi_v6/merged_<ARM>.bam(.bai)   real files
#          human/pm/02_gt_bams_umi_v6/T{1-4}_<ARM>.bam(.bai)   symlinks to the same per-teratoma BAMs
# Check  : merged read count == sum of the four per-teratoma counts (idxstats).
#
# Usage  : qsub -t 1-5 -q slice human/pm/scripts/merge_gt_umi_bams_v6.sh
# =============================================================================
#$ -N merge_umi_v6
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 4
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00

set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"

BASE="${SCRNA_DIR:-.}/human/pm"
SRC="$BASE/02_gt_bams_umi"
OUT="$BASE/02_gt_bams_umi_v6"
ARMS=(AAVS1_ctrl ADAR1_KO ADARB1_KO ADARB2_KO NTC_ctrl)
ARM="${ARMS[$((SGE_TASK_ID - 1))]}"
[ -n "$ARM" ] || { echo "ERROR: no arm for task ${SGE_TASK_ID}" >&2; exit 1; }
mkdir -p "$OUT"

reads() { samtools idxstats "$1" | awk '{s += $3 + $4} END {print s}'; }

IN=(); SUM=0
for T in T1 T2 T3 T4; do
  # PERSAMPLE_DIR (TSO-clipped rerun): read <dir>/<T>_<ARM>/02_dedup/<T>_<ARM>_umiDup.bam instead of $SRC
  if [ -n "${PERSAMPLE_DIR:-}" ]; then b=$(readlink -f "$PERSAMPLE_DIR/${T}_${ARM}/02_dedup/${T}_${ARM}_umiDup.bam"); else b=$(readlink -f "$SRC/${T}_${ARM}.bam"); fi
  [ -s "$b" ] && [ -s "${b}.bai" ] || { echo "ERROR: missing $b(.bai)" >&2; exit 1; }
  ln -sfn "$b" "$OUT/${T}_${ARM}.bam"
  ln -sfn "${b}.bai" "$OUT/${T}_${ARM}.bam.bai"
  IN+=("$b"); SUM=$((SUM + $(reads "$b")))
done

echo "=== $(date) === merge ${ARM}: ${#IN[@]} per-teratoma UMI BAMs, ${SUM} reads  ($(samtools --version | head -1)) ==="
TMP="$OUT/.merged_${ARM}.tmp.bam"
samtools merge -@ 4 -f -o "$TMP" "${IN[@]}"
samtools index -@ 4 "$TMP"
N=$(reads "$TMP")
[ "$N" -eq "$SUM" ] || { echo "ERROR: merged ${N} != sum ${SUM}" >&2; exit 1; }
mv -f "$TMP" "$OUT/merged_${ARM}.bam"
mv -f "${TMP}.bai" "$OUT/merged_${ARM}.bam.bai"
echo "=== $(date) === wrote merged_${ARM}.bam  (${N} reads = sum of T1-T4) ==="
