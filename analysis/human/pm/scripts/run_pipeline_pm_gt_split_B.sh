#!/bin/bash
# =============================================================================
# run_pipeline_pm_gt_split_B.sh  [sheet|run|collect]   (2026-10-13)
#
# Split-by-cells test: the v6 caller grid on NTC cell half B (split_gt_bams_halves.sh), all three tiers in
# one run: merged arm NTC_ctrl, T1-T4_NTC_ctrl, the 11 kept NTC cell types. Unit ids are the v6 ids so the
# scorers only need a different results dir. FRESH outdir / work-dir (results_gt_split_B, work_gt_split_B).
# Same parameters as run_pipeline_pm_gt_{merged_full,persample,celltype}_v6.sh, plus
#   --callers + cellsnp           (v6 ran cellSNP-lite separately, cellsnp_v6.sh; same CELLSNP_LITE + DOWNSTREAM_R)
#   --mismatch_window_filter true (the v6 tasks ran downstream_filter.R without the flag = its default 1;
#                                  downstream_r.nf now passes the param, default false)
#   --save_intermediates 4        (step-4 / step-5 tables in the work dirs = what prewindow_rerun_v6.sh recovered)
# Dedup: prebuilt sampleRmDup / sampleUmiDup = the half-B cut of the v6 per-GEM-well dedup BAMs.
#   bash run_pipeline_pm_gt_split_B.sh sheet
#   qsub -v STAGE=run run_pipeline_pm_gt_split_B.sh
#   qsub -v STAGE=collect run_pipeline_pm_gt_split_B.sh   (after the run)
# collect: <prefix>_step4_cov_edit_filtered / _step5_nomm.tsv.gz next to each published pre100 (empty file when
# the step had no sites), from the work dir whose pre100 is identical to the published one.
#$ -N nf_pm_gt_splitB
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 1
#$ -l h_vmem=64G
#$ -l h_rt=336:00:00
set -eo pipefail
export NXF_OPTS="-Xms1g -Xmx10g"
BASEDIR="${SCRNA_DIR:-.}"
NFDIR="${BASEDIR}/nextflow"
REFDIR="${BASEDIR}/human/ref_annots"
DATADIR="${BASEDIR}/human/pm"
S="$DATADIR/gt_eval_v8/split_cells"
SHEET="$S/samplesheet_split_B.csv"
OUTDIR="$DATADIR/results_gt_split_B"
WORK="$DATADIR/work_gt_split_B"
STAGE=${STAGE:-${1:-run}}

if [ "$STAGE" = "sheet" ]; then
  BC="$S/barcodes"; BM="$S/bams"
  row() {  # sample_id bam_stem barcode_file
    printf '%s,%s,%s,%s,10xv3,"","",NTC_ctrl,%s,%s\n' "$1" "$BM/raw/$2_B.bam" "$3" "$3" "$BM/rmdup/$2_B.bam" "$BM/umi/$2_B.bam"
  }
  {
    echo "sample_id,bam,celltype_annot,called_barcodes,platform,fastq_r1,fastq_r2,group,rmdup_bam,umidup_bam"
    row NTC_ctrl merged_NTC_ctrl "$BC/merged_NTC_ctrl_B.tsv"
    for n in 1 2 3 4; do row "T${n}_NTC_ctrl" "T${n}_NTC_ctrl" "$BC/T${n}_NTC_ctrl_B.tsv"; done
    while read -r u; do row "$u" "$u" "$BC/${u}_B.tsv"; done < "$DATADIR/gt_eval_v8/unit_cov/celltype_keep_min20.txt"
  } > "$SHEET"
  miss=$(awk -F, 'NR > 1 { for (i = 2; i <= 4; i++) if (system("test -s " $i)) print $i;
                           for (i = 9; i <= 10; i++) if (system("test -s " $i " && test -s " $i ".bai")) print $i }' "$SHEET")
  [ -z "$miss" ] || { echo "ERROR: missing:" >&2; echo "$miss" >&2; exit 1; }
  echo "$(( $(wc -l < "$SHEET") - 1 )) units -> $SHEET"
  exit 0
fi

if [ "$STAGE" = "run" ]; then
  RUNDIR="$DATADIR/nf_run_gt_split_B"; mkdir -p "$RUNDIR"; cd "$RUNDIR"
  export PS1="${PS1:-}"
  export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  source ~/.bashrc 2>/dev/null || true
  mamba activate scRNAEditing
  echo "=== $(date) === split-by-cells half B, full grid ==="
  nextflow run "${NFDIR}/main.nf" \
    -profile sge,conda \
    --input "$SHEET" \
    --outdir "$OUTDIR" \
    --genome_fasta "${REFDIR}/hg38.analysisSet.fa" \
    --genome_gtf "${REFDIR}/hg38.ncbiRefSeq.gtf" \
    --snp_file "${BASEDIR}/annotations/SNP_files/dbSNP156_common_GRanges.RData" \
    --strand_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --alu_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --rediportal "${REFDIR}/REDIportal_hg38_v2.vcf.gz" \
    --known_sites_vcf "${REFDIR}/dbSNP156_GRCh38.vcf.gz" \
    --refseq_bed "${BASEDIR}/review/ref_annots/hg38_RefSeq.bed" \
    --main_chroms "${REFDIR}/chrlist_main_chrs.txt" \
    --callers "bcftools,reditools3,gatk,jacusa2,scallele,cellsnp" \
    --scallele_conda "${CONDA_ENVS:-/opt/conda/envs}/scAllele" \
    --dedup "none,coordinate,umi" \
    --strand_modes "annotation,hybrid,caller" \
    --bam_types "sample,sampleRmDup,sampleUmiDup,recal,recalq20" \
    --qual 30 \
    --min_cov 10 \
    --requant_editing false \
    --percell_analyses false \
    --mismatch_window_filter true \
    --save_intermediates 4 \
    -work-dir "$WORK" \
    -resume
  echo "=== $(date) === Pipeline complete ==="
  exit 0
fi

if [ "$STAGE" = "collect" ]; then
  n_ok=0; n_empty=0; n_bad=0
  IDX=$(mktemp); trap 'rm -f "$IDX"' EXIT
  find "$WORK" -mindepth 3 -maxdepth 3 -name '*_pre100.tsv.gz' -type f > "$IDX"
  while read -r pub; do
    d=$(dirname "$pub"); prefix=$(basename "$pub" _pre100.tsv.gz)
    wd=""
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
  exit 0
fi
echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
