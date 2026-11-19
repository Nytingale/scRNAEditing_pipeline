#!/bin/bash
# aorta_coord_bams.sh STAGE=dedup|merge   (user 2026-11-12)
# Coordinate-deduplicated (MarkDuplicates) BAMs for the aorta Adar1KD mice, the knockdown units of a coordinate-built
# truth (cross-truth test; the control units are the pipeline's own results_v8/<unit>/02_dedup/<unit>_rmDup.bam).
# Mirrors the pipeline path of the control mice: aorta_pipeline_bams.sh sample (STARsolo BAM, Myh11+ cells by CB,
# @RG ID/SM = mouse) -> ADD_MD_TAG (calmd only if MD is missing) -> FILTER_UNIQUELY_MAPPED ([NH]==1) -> MARKDUP_CHROM
# (per main chromosome + _other, --BARCODE_TAG CR --REMOVE_DUPLICATES true; one mouse = one well) -> gather.
#   dedup  array 1-4: 1-3 Adar1KD_{1,2,3}_tdT -> pipeline_bams_coord/; 4 Control_3_tdT -> pipeline_bams_coord/check/,
#          compared per contig with results_v8/Control_3_tdT/02_dedup/Control_3_tdT_rmDup.bam (equivalence test)
#   merge  merged_KD.bam (3 read groups; read count == sum checked)
# Usage: j=$(qsub -terse -v STAGE=dedup -t 1-4 mouse/aorta_gse255661/scripts/aorta_coord_bams.sh | cut -d. -f1)
#        qsub -v STAGE=merge -hold_jid $j mouse/aorta_gse255661/scripts/aorta_coord_bams.sh
#$ -S /bin/bash
#$ -N aorta_coord
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/logs/
#$ -j y
#$ -pe threaded 8
#$ -l h_rt=24:00:00
set -euo pipefail
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
B=${SCRNA_DIR:-.}
AO=$B/mouse/aorta_gse255661
O=$AO/pipeline_bams_coord
FA=$B/mouse/ref_annots/mm39.fa
MAIN=$B/mouse/ref_annots/chrlist_main_chrs.txt
NT=${NSLOTS:-8}
STAGE=${STAGE:-${1:-}}
reads() { samtools idxstats "$1" | awk '{ s += $3 + $4 } END { print s }'; }

markdup_chunk() {  # <in.bam> <chunk> <workdir>
  local in=$1 c=$2 w=$3 regs
  if [ "$c" = _other ]; then
    regs=$(samtools idxstats "$in" | awk -F'\t' 'NR == FNR { m[$1] = 1; next } $1 != "*" && $3 > 0 && !($1 in m) { print $1 }' "$MAIN" - | tr '\n' ' ')
  else
    regs=$c
  fi
  if [ -z "$regs" ] || [ "$(samtools view -c "$in" $regs)" -eq 0 ]; then samtools view -H -b -o "$w/d_$c.bam" "$in"; return; fi
  samtools view -b -o "$w/t_$c.bam" "$in" $regs
  gatk --java-options "-Xmx5G" MarkDuplicates -I "$w/t_$c.bam" -O "$w/d_$c.bam" -M "$w/m_$c.txt" \
    --BARCODE_TAG CR --REMOVE_DUPLICATES true --QUIET true --VERBOSITY WARNING
  rm -f "$w/t_$c.bam"
}
export -f markdup_chunk
export MAIN

case "$STAGE" in
dedup)
  case $SGE_TASK_ID in
    1|2|3) S=Adar1KD_${SGE_TASK_ID}_tdT; OUT=$O ;;
    4) S=Control_3_tdT; OUT=$O/check ;;
  esac
  BAM=$AO/solo/$S/${S}_Aligned.sortedByCoord.out.bam; BC=$AO/truth/barcodes/$S.txt
  [ -s "$BAM" ] && [ -s "$BC" ] || { echo "ERROR: missing $BAM or $BC" >&2; exit 1; }
  W=$B/tmp/coord_bams/$S; rm -rf "$W"; mkdir -p "$W" "$OUT"
  echo "=== $(date) $S: $(wc -l < "$BC") Myh11+ cells"
  samtools view -@ "$NT" -u -D "CB:$BC" "$BAM" \
    | samtools addreplacerg -@ "$NT" -m overwrite_all -r "ID:$S" -r "SM:$S" -r "PL:ILLUMINA" -o "$W/sample.bam" -
  if [ "$(samtools view "$W/sample.bam" | head -1000 | grep -c 'MD:' || true)" -eq 0 ]; then
    echo "MD tag missing, adding with samtools calmd"
    samtools calmd -@ "$NT" -b "$W/sample.bam" "$FA" > "$W/md.bam" 2> /dev/null; mv "$W/md.bam" "$W/sample.bam"
  fi
  samtools view -b -e '[NH]==1' -@ "$NT" -o "$W/uniq.bam" "$W/sample.bam"
  samtools index -@ "$NT" "$W/uniq.bam"; rm -f "$W/sample.bam"
  { cat "$MAIN"; echo _other; } > "$W/chunks.txt"
  xargs -P "$NT" -I{} bash -c 'markdup_chunk "$0" "$1" "$2"' "$W/uniq.bam" {} "$W" < "$W/chunks.txt"
  for c in $(cat "$W/chunks.txt"); do [ -e "$W/d_$c.bam" ] || { echo "ERROR: missing chunk $c" >&2; exit 1; }; done
  samtools merge -@ "$NT" -c -p -f -o "$OUT/$S.bam" $(sed "s|^|$W/d_|; s|$|.bam|" "$W/chunks.txt")
  samtools index -@ "$NT" "$OUT/$S.bam"
  cat "$W"/m_*.txt > "$OUT/$S.dup_metrics.txt" 2> /dev/null || : > "$OUT/$S.dup_metrics.txt"
  printf "%s\t%s\t%s\n" "$S" "$(reads "$W/uniq.bam")" "$(reads "$OUT/$S.bam")" > "$OUT/$S.counts.tsv"
  echo "=== $(date) $S: $(cut -f2,3 "$OUT/$S.counts.tsv") (unique-mapped -> kept reads)"
  if [ "$SGE_TASK_ID" = 4 ]; then
    REF=$AO/results_v8/$S/02_dedup/${S}_rmDup.bam
    if diff <(samtools idxstats "$OUT/$S.bam" | cut -f1,3) <(samtools idxstats "$REF" | cut -f1,3) > "$OUT/$S.idxstats_diff.txt"; then
      echo "EQUIVALENT: per-contig read counts identical to $REF"
    else
      echo "NOT EQUIVALENT: $(wc -l < "$OUT/$S.idxstats_diff.txt") diff lines vs $REF" >&2; exit 1
    fi
  fi
  rm -rf "$W"
  ;;
merge)
  IN="$O/Adar1KD_1_tdT.bam $O/Adar1KD_2_tdT.bam $O/Adar1KD_3_tdT.bam"
  [ -s "$O/check/Control_3_tdT.idxstats_diff.txt" ] && { echo "ERROR: Control_3 equivalence test failed" >&2; exit 1; }
  [ -e "$O/check/Control_3_tdT.bam" ] || { echo "ERROR: Control_3 equivalence test did not run" >&2; exit 1; }
  samtools merge -@ "$NT" -c -p -f -o "$O/merged_KD.bam" $IN
  samtools index -@ "$NT" "$O/merged_KD.bam"
  s=0; for f in $IN; do s=$(( s + $(reads "$f") )); done
  m=$(reads "$O/merged_KD.bam")
  [ "$s" = "$m" ] || { echo "ERROR: merged $m != sum $s" >&2; exit 1; }
  echo "=== $(date) merged_KD.bam: $m reads (= sum of inputs), $(samtools view -H "$O/merged_KD.bam" | grep -c '^@RG') read groups"
  ;;
*) echo "STAGE=dedup|merge" >&2; exit 1 ;;
esac
