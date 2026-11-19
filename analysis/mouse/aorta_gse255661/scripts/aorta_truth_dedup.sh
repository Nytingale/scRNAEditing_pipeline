#!/bin/bash
# =============================================================================
# aorta_truth_dedup.sh   STAGE=dedup (SGE array 1-6, one task per mouse) | STAGE=merge (array 1-2, one per arm)
#
# Ground-truth BAMs for GSE255661, 3 Adar1KD vs 3 Control tdT+ mice (10 weeks):
#   dedup: STARsolo BAM -> Myh11+ cells only (truth/barcodes/<sample>.txt, CB tag)
#          -> reads with an uncorrectable UMI (UB:Z:-) dropped -> umi_tools dedup
#          per chromosome with the pipeline settings (--per-cell --umi-tag=UB
#          --method=unique; cell tag CB, the tag the cell list is defined on).
#          One mouse = one GEM well, so each sample is its own dedup unit.
#   merge: per-mouse UMI BAMs -> merged_<ARM>.bam (read count == sum checked).
#
# Inputs : solo/<sample>/<sample>_Aligned.sortedByCoord.out.bam, truth/design.tsv
# Outputs: truth/bams/<sample>.bam(.bai), truth/bams/merged_{KD,CTRL}.bam(.bai)
# Usage  : STAGE=dedup qsub -t 1-6 scripts/aorta_truth_dedup.sh
#          STAGE=merge qsub -t 1-2 -hold_jid <dedup job> scripts/aorta_truth_dedup.sh
# =============================================================================
#$ -N aorta_dedup
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/logs/
#$ -j y
#$ -V
#$ -pe threaded 8
#$ -l mem=8G
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00

set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:/usr/local/bin:/usr/bin:/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"

BASE="${SCRNA_DIR:-.}/mouse/aorta_gse255661"
OUT="${OUT:-$BASE/truth/bams}"   # OUT=$BASE/truth/bams_seed8787: seeded rebuild (2026-11-11)
STAGE="${STAGE:?STAGE=dedup|merge}"
NCPU="${NSLOTS:-8}"
mkdir -p "$OUT"

reads() { samtools idxstats "$1" | awk '{s += $3 + $4} END {print s}'; }

if [ "$STAGE" = "dedup" ]; then
  S=$(awk -F'\t' -v i="$SGE_TASK_ID" 'NR == i + 1 { print $1 }' "$BASE/truth/design.tsv")
  [ -n "$S" ] || { echo "ERROR: no sample for task $SGE_TASK_ID" >&2; exit 1; }
  BAM="$BASE/solo/$S/${S}_Aligned.sortedByCoord.out.bam"
  BC="$BASE/truth/barcodes/$S.txt"
  [ -s "$BAM" ] && [ -s "$BC" ] || { echo "ERROR: missing $BAM or $BC" >&2; exit 1; }
  W="${SCRNA_DIR:-.}/tmp/aorta_dedup_$(basename "$OUT")/$S"
  rm -rf "$W"; mkdir -p "$W"
  echo "=== $(date) === $S: $(wc -l < "$BC") Myh11+ cells  ($(samtools --version | head -1); umi_tools $(umi_tools --version 2>&1 | awk '{print $NF}'))"

  samtools idxstats "$BAM" | awk '$3 > 0 && $1 ~ /^chr([0-9]+|X|Y|M)$/ { print $1 }' > "$W/chroms.txt"
  dedup_chrom() {
    local c=$1
    samtools view -b -D "CB:$BC" -e '[UB] != "-"' -o "$W/$c.in.bam" "$BAM" "$c"
    samtools index "$W/$c.in.bam"
    umi_tools dedup --per-cell --cell-tag=CB --umi-tag=UB --extract-umi-method=tag --method=unique --random-seed=${SEED:-8787} \
      -I "$W/$c.in.bam" -S "$W/$c.dd.bam" --log="$W/$c.log"
    echo "$c $(samtools view -c "$W/$c.in.bam") -> $(samtools view -c "$W/$c.dd.bam")"
    rm -f "$W/$c.in.bam" "$W/$c.in.bam.bai"
  }
  export -f dedup_chrom; export BAM BC W
  xargs -P "$NCPU" -I{} bash -c 'set -eo pipefail; dedup_chrom {}' < "$W/chroms.txt" | sort -V > "$W/counts.txt"
  [ "$(wc -l < "$W/counts.txt")" -eq "$(wc -l < "$W/chroms.txt")" ] || { echo "ERROR: not every chromosome finished" >&2; exit 1; }
  cat "$W/counts.txt"

  samtools merge -@ "$NCPU" -c -p -f -o "$OUT/.$S.tmp.bam" $(sed "s#.*#$W/&.dd.bam#" "$W/chroms.txt")
  samtools index -@ "$NCPU" "$OUT/.$S.tmp.bam"
  N=$(reads "$OUT/.$S.tmp.bam"); E=$(awk '{s += $4} END {print s}' "$W/counts.txt")
  [ "$N" -eq "$E" ] || { echo "ERROR: merged $N != per-chromosome sum $E" >&2; exit 1; }
  mv -f "$OUT/.$S.tmp.bam" "$OUT/$S.bam"; mv -f "$OUT/.$S.tmp.bam.bai" "$OUT/$S.bam.bai"
  cat "$W"/*.log > "$OUT/$S.umi_dedup.log"
  awk -v s="$S" '{ i += $2; o += $4 } END { printf "%s\t%d\t%d\t%.3f\n", s, i, o, o / i }' "$W/counts.txt" > "$OUT/$S.counts.tsv"
  rm -rf "$W"
  echo "=== $(date) === wrote $OUT/$S.bam ($N reads)"

elif [ "$STAGE" = "merge" ]; then
  ARMS=(KD CTRL); ARM="${ARMS[$((SGE_TASK_ID - 1))]}"
  [ -n "$ARM" ] || { echo "ERROR: no arm for task $SGE_TASK_ID" >&2; exit 1; }
  IN=(); SUM=0
  for S in $(awk -F'\t' -v a="$ARM" 'NR > 1 && $2 == a { print $1 }' "$BASE/truth/design.tsv"); do
    b="$OUT/$S.bam"; [ -s "$b" ] && [ -s "$b.bai" ] || { echo "ERROR: missing $b(.bai)" >&2; exit 1; }
    IN+=("$b"); SUM=$((SUM + $(reads "$b")))
  done
  [ "${#IN[@]}" -eq 3 ] || { echo "ERROR: expected 3 mice in $ARM, found ${#IN[@]}" >&2; exit 1; }
  echo "=== $(date) === merge $ARM: ${IN[*]} ($SUM reads)"
  samtools merge -@ "$NCPU" -f -o "$OUT/.merged_$ARM.tmp.bam" "${IN[@]}"
  samtools index -@ "$NCPU" "$OUT/.merged_$ARM.tmp.bam"
  N=$(reads "$OUT/.merged_$ARM.tmp.bam")
  [ "$N" -eq "$SUM" ] || { echo "ERROR: merged $N != sum $SUM" >&2; exit 1; }
  mv -f "$OUT/.merged_$ARM.tmp.bam" "$OUT/merged_$ARM.bam"; mv -f "$OUT/.merged_$ARM.tmp.bam.bai" "$OUT/merged_$ARM.bam.bai"
  echo "=== $(date) === wrote merged_$ARM.bam ($N reads)"
else
  echo "ERROR: STAGE must be dedup or merge" >&2; exit 1
fi
