#!/bin/bash
# umimd_bams.sh STAGE=dedup|merge|splitprep
# UMI-aware MarkDuplicates test: same grouping as the UMI-dedup BAMs (cell, UMI, 5' position, strand; reads with UB:Z:-
# dropped) but the kept read is MarkDuplicates' highest base-quality sum instead of umi_tools' random pick.
# Tag XU = <cell tag><UB> (DNA only; MarkDuplicates --BARCODE_TAG rejects other characters). Cell tag as the existing
# UMI BAMs: teratoma CR (pipeline umi_tools --cell-tag CR), aorta CB (aorta_truth_dedup.sh --cell-tag=CB).
#   dedup      array 1-7: 1-4 teratoma T{1-4}_NTC_ctrl (02_gt_bams/), 5-7 aorta control mice (pipeline_bams/);
#              8-15 teratoma T{1-4}_ADAR1_KO, T{1-4}_ADARB1_KO and 16-18 aorta Adar1KD mice (STARsolo BAM, Myh11+ cells as
#              aorta_pipeline_bams.sh): the knockout units of a umimd-built truth;
#              per chromosome in parallel, then concatenated in header order
#   merge      array 1-5: 1 teratoma merged_NTC_ctrl (T1-T4), 2 aorta merged_CTRL (3 mice), 3 merged_ADAR1_KO,
#              4 merged_ADARB1_KO, 5 aorta merged_KD
#   splitprep  login node: plan dir for split_gt_bams_celltype.sh (NTC unit, new per-teratoma BAMs)
# Outputs: human/pm/02_gt_bams_umimd/, mouse/aorta_gse255661/pipeline_bams_umimd/
#$ -S /bin/bash
#$ -N umimd
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 8
#$ -l h_rt=24:00:00
set -euo pipefail
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
B=${SCRNA_DIR:-.}
PM=$B/human/pm
AO=$B/mouse/aorta_gse255661
OPM=$PM/02_gt_bams_umimd
OAO=$AO/pipeline_bams_umimd
STAGE=${STAGE:-${1:-}}
NT=${NSLOTS:-8}

dedup_chrom() {  # <in.bam> <chrom> <cell tag> <workdir> [cell barcode list]
  local in=$1 c=$2 ct=$3 w=$4 bc=${5:-}
  samtools view -h ${bc:+-D "$ct:$bc"} "$in" "$c" | awk -v ct="$ct:Z:" 'BEGIN { OFS = "\t" } /^@/ { print; next }
    { cb = ""; ub = ""; for (i = 12; i <= NF; i++) { t = substr($i, 1, 5); if (t == "UB:Z:") ub = substr($i, 6); else if (t == ct) cb = substr($i, 6) }
      if (ub == "" || ub == "-" || cb == "" || cb == "-") next
      print $0, "XU:Z:" cb ub }' | samtools view -b -o "$w/t_$c.bam" -
  gatk --java-options "-Xmx5G" MarkDuplicates -I "$w/t_$c.bam" -O "$w/d_$c.bam" -M "$w/m_$c.txt" \
    --BARCODE_TAG XU --REMOVE_DUPLICATES true --QUIET true --VERBOSITY WARNING
  rm -f "$w/t_$c.bam"
}
export -f dedup_chrom

case "$STAGE" in
dedup)
  case $SGE_TASK_ID in
    1|2|3|4) S=T${SGE_TASK_ID}_NTC_ctrl; IN=$PM/02_gt_bams/$S.bam; CT=CR; O=$OPM ;;
    5) S=Control_1_tdT_nontdT_combined; IN=$AO/pipeline_bams/$S.bam; CT=CB; O=$OAO ;;
    6) S=Control_2_tdT_nontdT_combined; IN=$AO/pipeline_bams/$S.bam; CT=CB; O=$OAO ;;
    7) S=Control_3_tdT; IN=$AO/pipeline_bams/$S.bam; CT=CB; O=$OAO ;;
    8|9|10|11) S=T$((SGE_TASK_ID - 7))_ADAR1_KO; IN=$PM/02_gt_bams/$S.bam; CT=CR; O=$OPM ;;
    12|13|14|15) S=T$((SGE_TASK_ID - 11))_ADARB1_KO; IN=$PM/02_gt_bams/$S.bam; CT=CR; O=$OPM ;;
    16|17|18) S=Adar1KD_$((SGE_TASK_ID - 15))_tdT; IN=$AO/solo/$S/${S}_Aligned.sortedByCoord.out.bam; BC=$AO/truth/barcodes/$S.txt; CT=CB; O=$OAO ;;
  esac
  W=$B/tmp/umimd/$S; rm -rf "$W"; mkdir -p "$W" "$O"
  samtools idxstats "$IN" | awk '$1 != "*" && $3 > 0 { print $1 }' > "$W/chroms.txt"
  echo "=== $(date) $S: $(wc -l < "$W/chroms.txt") contigs, cell tag $CT"
  xargs -P "$NT" -I{} bash -c 'dedup_chrom "$0" "$1" "$2" "$3" "$4"' "$IN" {} "$CT" "$W" "${BC:-}" < "$W/chroms.txt"
  for c in $(cat "$W/chroms.txt"); do [ -s "$W/d_$c.bam" ] || { echo "ERROR: missing $c" >&2; exit 1; }; done
  samtools cat -o "$O/$S.bam" $(sed "s|^|$W/d_|; s|$|.bam|" "$W/chroms.txt")
  samtools index "$O/$S.bam"
  cat "$W"/m_*.txt > "$O/$S.dup_metrics.txt"
  n_in=$(samtools view -c -F 0x904 ${BC:+-D "CB:$BC"} "$IN"); n_out=$(samtools view -c -F 0x904 "$O/$S.bam")
  printf "%s\t%s\t%s\n" "$S" "$n_in" "$n_out" > "$O/$S.counts.tsv"
  echo "=== $(date) $S done: $n_in -> $n_out primary reads"
  rm -rf "$W"
  ;;
merge)
  case $SGE_TASK_ID in
    1) O=$OPM; M=$O/merged_NTC_ctrl.bam; IN="$O/T1_NTC_ctrl.bam $O/T2_NTC_ctrl.bam $O/T3_NTC_ctrl.bam $O/T4_NTC_ctrl.bam" ;;
    2) O=$OAO; M=$O/merged_CTRL.bam; IN="$O/Control_1_tdT_nontdT_combined.bam $O/Control_2_tdT_nontdT_combined.bam $O/Control_3_tdT.bam" ;;
    3|4) O=$OPM; a=ADAR1_KO; [ "$SGE_TASK_ID" = 4 ] && a=ADARB1_KO; M=$O/merged_$a.bam; IN="$O/T1_$a.bam $O/T2_$a.bam $O/T3_$a.bam $O/T4_$a.bam" ;;
    5) O=$OAO; M=$O/merged_KD.bam; IN="$O/Adar1KD_1_tdT.bam $O/Adar1KD_2_tdT.bam $O/Adar1KD_3_tdT.bam" ;;
  esac
  samtools merge -@ "$NT" -c -p -f -o "$M" $IN
  samtools index "$M"
  s=0; for f in $IN; do s=$(( s + $(samtools idxstats "$f" | awk '{ t += $3 + $4 } END { print t }') )); done
  m=$(samtools idxstats "$M" | awk '{ t += $3 + $4 } END { print t }')
  [ "$s" = "$m" ] || { echo "ERROR: merged $m != sum $s" >&2; exit 1; }
  echo "=== $(date) $M: $m reads (= sum of inputs)"
  ;;
splitprep)
  P=$PM/_celltype_frags_umimd; rm -rf "$P"; mkdir -p "$P"
  S6=$PM/_celltype_frags_umi_v6
  cp "$S6/_barcode_prefix_map.tsv" "$S6/_cells_NTC_ctrl.tsv" "$S6/_keep_NTC_ctrl.txt" "$P/"
  sed "s|/02_gt_bams_umi_v6/|/02_gt_bams_umimd/|" "$S6/_plan_NTC_ctrl.tsv" > "$P/_plan_NTC_ctrl.tsv"
  echo NTC_ctrl > "$P/_expected.txt"
  [ -d "$S6/barcodes" ] && cp -r "$S6/barcodes" "$P/"
  cat "$P/_plan_NTC_ctrl.tsv"
  ;;
*) echo "STAGE=dedup|merge|splitprep" >&2; exit 1 ;;
esac
