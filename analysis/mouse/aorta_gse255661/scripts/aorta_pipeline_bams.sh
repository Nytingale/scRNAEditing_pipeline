#!/bin/bash
# =============================================================================
# aorta_pipeline_bams.sh   STAGE=sample (SGE array 1-3, one per control mouse) | STAGE=merge
#
# Caller-input BAMs for the aorta benchmark (controls only; the KD arm is only used by the truth):
#   sample: STARsolo BAM -> reads of the Myh11+ cells (truth/barcodes/<sample>.txt, CB tag), no dedup
#           -> @RG ID/SM = <sample> on every read, so the pipeline's coordinate / UMI dedup runs per mouse
#           (one mouse = one GEM well) on the merged BAM.
#   merge : the 3 control mice -> merged_CTRL.bam (read count == sum checked), 3 read groups.
# UMI-dedup BAMs are the truth BAMs (truth/bams/, aorta_truth_dedup.sh), passed as umidup_bam, as the
# teratoma used 02_gt_bams_umi_v6 for both truth and callers.
# Also writes the samplesheet pipeline_bams/samplesheet_aorta_ctrl.csv (4 rows: CTRL arm + 3 mice).
# Usage: STAGE=sample qsub -t 1-3 scripts/aorta_pipeline_bams.sh
#        STAGE=merge  qsub -hold_jid <sample job> scripts/aorta_pipeline_bams.sh
# =============================================================================
#$ -N aorta_pbam
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/logs/
#$ -j y
#$ -V
#$ -pe threaded 4
#$ -l h_vmem=4G
#$ -l h_rt=24:00:00
set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:/usr/local/bin:/usr/bin:/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"

BASE="${SCRNA_DIR:-.}/mouse/aorta_gse255661"
OUT="$BASE/pipeline_bams"
STAGE="${STAGE:?STAGE=sample|merge}"
NCPU="${NSLOTS:-4}"
mkdir -p "$OUT"
reads() { samtools idxstats "$1" | awk '{s += $3 + $4} END {print s}'; }
CTRL=($(awk -F'\t' 'NR > 1 && $2 == "CTRL" { print $1 }' "$BASE/truth/design.tsv"))
[ "${#CTRL[@]}" -eq 3 ] || { echo "ERROR: expected 3 control mice" >&2; exit 1; }

if [ "$STAGE" = "sample" ]; then
  S="${CTRL[$((SGE_TASK_ID - 1))]}"
  BAM="$BASE/solo/$S/${S}_Aligned.sortedByCoord.out.bam"; BC="$BASE/truth/barcodes/$S.txt"
  [ -s "$BAM" ] && [ -s "$BC" ] || { echo "ERROR: missing $BAM or $BC" >&2; exit 1; }
  echo "=== $(date) === $S: $(wc -l < "$BC") Myh11+ cells"
  samtools view -@ "$NCPU" -u -D "CB:$BC" "$BAM" \
    | samtools addreplacerg -@ "$NCPU" -m overwrite_all -r "ID:$S" -r "SM:$S" -r "PL:ILLUMINA" -o "$OUT/.$S.tmp.bam" -
  samtools index -@ "$NCPU" "$OUT/.$S.tmp.bam"
  mv -f "$OUT/.$S.tmp.bam" "$OUT/$S.bam"; mv -f "$OUT/.$S.tmp.bam.bai" "$OUT/$S.bam.bai"
  awk -v OFS='\t' '{ print $1, "Myh11_pos" }' "$BC" > "$OUT/$S.barcodes.tsv"
  echo "=== $(date) === wrote $S.bam ($(reads "$OUT/$S.bam") reads; @RG: $(samtools view -H "$OUT/$S.bam" | grep -c '^@RG'))"

elif [ "$STAGE" = "merge" ]; then
  IN=(); SUM=0
  for S in "${CTRL[@]}"; do
    b="$OUT/$S.bam"; [ -s "$b" ] && [ -s "$b.bai" ] || { echo "ERROR: missing $b" >&2; exit 1; }
    IN+=("$b"); SUM=$((SUM + $(reads "$b")))
  done
  samtools merge -@ "$NCPU" -f -o "$OUT/.merged_CTRL.tmp.bam" "${IN[@]}"
  samtools index -@ "$NCPU" "$OUT/.merged_CTRL.tmp.bam"
  N=$(reads "$OUT/.merged_CTRL.tmp.bam")
  [ "$N" -eq "$SUM" ] || { echo "ERROR: merged $N != sum $SUM" >&2; exit 1; }
  NRG=$(samtools view -H "$OUT/.merged_CTRL.tmp.bam" | grep -c '^@RG')
  [ "$NRG" -eq 3 ] || { echo "ERROR: merged BAM has $NRG read groups, expected 3" >&2; exit 1; }
  mv -f "$OUT/.merged_CTRL.tmp.bam" "$OUT/merged_CTRL.bam"; mv -f "$OUT/.merged_CTRL.tmp.bam.bai" "$OUT/merged_CTRL.bam.bai"
  cat $(printf "$OUT/%s.barcodes.tsv " "${CTRL[@]}") > "$OUT/merged_CTRL.barcodes.tsv"
  {
    echo "sample_id,bam,celltype_annot,called_barcodes,platform,fastq_r1,fastq_r2,umidup_bam"
    echo "CTRL,$OUT/merged_CTRL.bam,$OUT/merged_CTRL.barcodes.tsv,$OUT/merged_CTRL.barcodes.tsv,10xv3,\"\",\"\",$BASE/truth/bams/merged_CTRL.bam"
    for S in "${CTRL[@]}"; do
      echo "$S,$OUT/$S.bam,$OUT/$S.barcodes.tsv,$OUT/$S.barcodes.tsv,10xv3,\"\",\"\",$BASE/truth/bams/$S.bam"
    done
  } > "$OUT/samplesheet_aorta_ctrl.csv"
  echo "=== $(date) === wrote merged_CTRL.bam ($N reads, $NRG read groups) and samplesheet_aorta_ctrl.csv"
else
  echo "ERROR: STAGE must be sample or merge" >&2; exit 1
fi
