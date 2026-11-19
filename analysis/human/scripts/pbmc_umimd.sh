#!/bin/bash
# pbmc_umimd.sh DS=jzlxv3|jzlcel STAGE=dedup|sheet|run|collect
# UMI-aware MarkDuplicates on the Ding PBMC UMI datasets (10x 3' v3, CEL-seq2), for the truth-free comparison with
# umi_tools and coordinate dedup. Same read grouping as the published umi_tools BAMs (results_v8/<s>/02_dedup/*_umiDup.bam:
# pipeline UMI_DEDUP_CHROM on results_v8/<s>/01_merged/<s>_merged.bam, --cell-tag CR, reads with UB:Z:- removed; one GEM
# well per sample) but the kept read is MarkDuplicates' highest base-quality sum (as human/pm/scripts/umimd_bams.sh).
#   dedup    array 1-N over the samples of DS: per chromosome in parallel, concatenated in header order
#            -> human/<DS>/02_bams_umimd/<s>.bam (+ .dup_metrics.txt, .counts.tsv)
#   sheet    login node: samplesheet_<DS>.csv rows + umidup_bam column -> samplesheet_<DS>_umimd.csv
#   run      published results_v8 settings restricted to --dedup umi --bam_types sampleUmiDup -> results_umimd/
#   collect  step-4 / step-5 tables next to each pre100 (as umimd_run.sh STAGE=collect)
# Usage: qsub -v DS=jzlxv3,STAGE=dedup -t 1-1 human/scripts/pbmc_umimd.sh ; DS=jzlxv3 STAGE=sheet bash ... ;
#        qsub -v DS=jzlxv3,STAGE=run human/scripts/pbmc_umimd.sh (submits collect when the pipeline ends)
#$ -S /bin/bash
#$ -N pbmc_umimd
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 8
#$ -l h_rt=48:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
DS=${DS:?DS=jzlxv3|jzlcel}
STAGE=${STAGE:?STAGE=dedup|sheet|run|collect}
D=$B/human/$DS
SRC_SHEET=$D/samplesheet_$DS.csv
SHEET=$D/samplesheet_${DS}_umimd.csv
O=$D/02_bams_umimd
OUTDIR=$D/results_umimd
WORK=$D/work_umimd
RUNDIR=$D/nf_run_umimd
SAMPLES=($(tail -n +2 "$SRC_SHEET" | cut -d, -f1 | awk '!s[$0]++'))
SELF=$B/human/scripts/pbmc_umimd.sh
QSUB=$(command -v qsub || echo /cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub)

dedup_chrom() {  # <in.bam> <chrom> <cell tag> <workdir>
  local in=$1 c=$2 ct=$3 w=$4
  samtools view -h "$in" "$c" | awk -v ct="$ct:Z:" 'BEGIN { OFS = "\t" } /^@/ { print; next }
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
  export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
  NT=${NSLOTS:-8}
  S=${SAMPLES[$(( SGE_TASK_ID - 1 ))]}
  IN=$D/results_v8/$S/01_merged/${S}_merged.bam
  [ -s "$IN" ] && [ -s "$IN.bai" ] || { echo "missing $IN" >&2; exit 1; }
  W=$B/tmp/umimd/${DS}_$S; rm -rf "$W"; mkdir -p "$W" "$O"
  n_sm=$(samtools view -H "$IN" | awk -F'\t' '/^@RG/ { for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) s[substr($i, 4)] = 1 } END { print length(s) }')
  [ "$n_sm" -le 1 ] || { echo "ERROR: $IN pools $n_sm GEM wells; per-well dedup not implemented here" >&2; exit 1; }
  samtools idxstats "$IN" | awk '$1 != "*" && $3 > 0 { print $1 }' > "$W/chroms.txt"
  echo "=== $(date) $DS $S: $(wc -l < "$W/chroms.txt") contigs, cell tag CR"
  xargs -P "$NT" -I{} bash -c 'dedup_chrom "$0" "$1" "$2" "$3"' "$IN" {} CR "$W" < "$W/chroms.txt"
  for c in $(cat "$W/chroms.txt"); do [ -s "$W/d_$c.bam" ] || { echo "ERROR: missing $c" >&2; exit 1; }; done
  samtools cat -o "$O/$S.bam" $(sed "s|^|$W/d_|; s|$|.bam|" "$W/chroms.txt")
  samtools index "$O/$S.bam"
  cat "$W"/m_*.txt > "$O/$S.dup_metrics.txt"
  P=$D/results_v8/$S/02_dedup/${S}_umiDup.bam
  n_in=$(samtools view -c -F 0x904 "$IN"); n_out=$(samtools view -c -F 0x904 "$O/$S.bam"); n_umi=$(samtools view -c -F 0x904 "$P")
  printf "sample\tinput\tumimd\tumi_tools\n%s\t%s\t%s\t%s\n" "$S" "$n_in" "$n_out" "$n_umi" > "$O/$S.counts.tsv"
  echo "=== $(date) $S done: $n_in -> $n_out primary reads (umi_tools: $n_umi)"
  rm -rf "$W"
  ;;
sheet)
  awk -F, -v OFS=, -v o="$O" 'NR == 1 { print $0, "umidup_bam"; next } { print $0, o "/" $1 ".bam" }' "$SRC_SHEET" > "$SHEET"
  echo "$SHEET: $(tail -n +2 "$SHEET" | wc -l) rows"
  for s in "${SAMPLES[@]}"; do [ -s "$O/$s.bam" ] && [ -s "$O/$s.bam.bai" ] && echo "  $s ok" || echo "  $s pending"; done
  ;;
run)
  [ -s "$SHEET" ] || { echo "missing $SHEET (STAGE=sheet)" >&2; exit 1; }
  for s in "${SAMPLES[@]}"; do [ -s "$O/$s.bam" ] && [ -s "$O/$s.bam.bai" ] || { echo "missing BAM $O/$s.bam" >&2; exit 1; }; done
  export NXF_OPTS="-Xms1g -Xmx10g"
  mkdir -p "$RUNDIR"; cd "$RUNDIR"
  export PS1="${PS1:-}"
  export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  set +u; source ~/.bashrc 2>/dev/null || true; mamba activate scRNAEditing; set -u
  R=$B/human/ref_annots
  echo "=== $(date) === umimd caller run, $DS -> $OUTDIR ==="
  nextflow run "$B/nextflow/main.nf" -profile sge,conda --input "$SHEET" --outdir "$OUTDIR" \
    --genome_fasta "$R/hg38.analysisSet.fa" --genome_gtf "$R/hg38.ncbiRefSeq.gtf" \
    --snp_file "$B/annotations/SNP_files/dbSNP156_common_GRanges.RData" \
    --strand_annots "$B/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --alu_annots "$B/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --rediportal "$R/REDIportal_hg38_v2.vcf.gz" --known_sites_vcf "$R/dbSNP156_GRCh38.vcf.gz" \
    --refseq_bed "$B/review/ref_annots/hg38_RefSeq.bed" --main_chroms "$R/chrlist_main_chrs.txt" \
    --callers "bcftools,reditools3,gatk,jacusa2,scallele,cellsnp" --scallele_conda "${CONDA_ENVS:-/opt/conda/envs}/scAllele" \
    --dedup "umi" --strand_modes "annotation,hybrid,caller" --bam_types "sampleUmiDup" \
    --qual 30 --min_cov 10 --bcftools_baq true --requant_editing false --percell_analyses false \
    --mismatch_window_filter true --clustered_filter false --save_intermediates 4 \
    -work-dir "$WORK" -resume
  echo "=== $(date) === Pipeline complete ==="
  $QSUB -v DS="$DS",STAGE=collect -pe threaded 1 "$SELF"
  ;;
collect)
  n_ok=0; n_empty=0; n_bad=0
  IDX=$(mktemp -p "$B/tmp"); trap 'rm -f "$IDX"' EXIT
  find "$WORK" -mindepth 3 -maxdepth 3 -name '*_pre100.tsv.gz' -type f > "$IDX"
  while read -r pub; do
    d=$(dirname "$pub"); prefix=$(basename "$pub" _pre100.tsv.gz); wd=""
    for c in $(grep -F "/${prefix}_pre100.tsv.gz" "$IDX" | xargs -r -n1 dirname); do
      cmp -s "$c/${prefix}_pre100.tsv.gz" "$pub" && { wd=$c; break; }
    done
    [ -n "$wd" ] || { echo "NO WORK DIR: $pub" >&2; n_bad=$((n_bad + 1)); continue; }
    for st in "step4_cov_edit_filtered:step4_cov_edit_filtered" "step5_variable_regions_removed:step5_nomm"; do
      src="$wd/${prefix}_filtered_${st%%:*}.tsv.gz"; dst="$d/${prefix}_${st##*:}.tsv.gz"
      if [ -s "$src" ]; then cp "$src" "$dst"; else : > "$dst"; n_empty=$((n_empty + 1)); fi
    done
    n_ok=$((n_ok + 1))
  done < <(find "$OUTDIR" -path '*/04_filtered/*' -name '*_pre100.tsv.gz')
  echo "collected $n_ok call sets ($n_empty empty step tables); $n_bad without a matching work dir"
  [ "$n_bad" -eq 0 ]
  ;;
*) echo "STAGE=dedup|sheet|run|collect" >&2; exit 1 ;;
esac
