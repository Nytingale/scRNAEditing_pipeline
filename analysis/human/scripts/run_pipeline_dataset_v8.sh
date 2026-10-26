#!/bin/bash
# =============================================================================
# run_pipeline_dataset_v8.sh  DATASET=<name>  STAGE=run|collect   (2026-10-18)
#
# Full configuration grid on a dataset with the teratoma v8 settings (human/pm/scripts/run_pipeline_pm_gt_v8.sh),
# for ranking configurations by the A-to-G composite rank (user 2026-10-18: every dataset, now and future, goes
# through all configurations). Add a dataset to the case block below.
#   jzlxv3 : Ding et al. PBMC1, 10x 3' v3 (4 libraries)              UMI platform
#   jzlss2 : Ding et al. PBMC1 + PBMC2 Smart-seq2, all 384 + 384 sequenced cells (per-cell BAMs; sample = all cells,
#            goodCells via qc_barcodes = the authors' 253 + 273 QC-passing cells; no UMIs: no umi dedup / *UmiDup)
# Datasets other than the ground truths (teratoma, Replogle) get every bam type: the sample-level tiers plus goodCells,
# annotated and celltype (each unprocessed, RmDup and, on UMI platforms, UmiDup). Cell-type call sets are published
# under <sample>__<celltype>/.
# FRESH outdir / work-dir / launch dir: <datadir>/{results_v8, work_v8, nf_run_v8}.
# Do not install into the scRNAEditing env while these run (invalidates -resume).
#   qsub -v DATASET=jzlxv3 human/scripts/run_pipeline_dataset_v8.sh
#   qsub -v DATASET=jzlxv3,STAGE=collect human/scripts/run_pipeline_dataset_v8.sh   (after the run)
# =============================================================================
#$ -N nf_ds_v8
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 1
#$ -l h_vmem=64G
#$ -l h_rt=336:00:00
set -euo pipefail
export NXF_OPTS="-Xms1g -Xmx10g"
BASEDIR="${SCRNA_DIR:-.}"
NFDIR="${BASEDIR}/nextflow"
REFDIR="${BASEDIR}/human/ref_annots"
DATASET=${DATASET:?DATASET=jzlxv3|jzlss2}
STAGE=${STAGE:-run}

case "$DATASET" in
  jzlxv3) DATADIR="$BASEDIR/human/jzlxv3"; SHEET="$DATADIR/samplesheet_jzlxv3.csv"; UMI=1 ;;
  jzlss2) DATADIR="$BASEDIR/human/jzlss2"; SHEET="$DATADIR/samplesheet_jzlss2_all768.csv"; UMI=0 ;;
  *) echo "ERROR: unknown DATASET '$DATASET'" >&2; exit 1 ;;
esac
BAMTYPES="sample,sampleRmDup,recal,recalq20,goodCells,goodCellsRmDup,annotated,annotatedRmDup,celltype,celltypeRmDup"
DEDUP="none,coordinate"
if [ "$UMI" -eq 1 ]; then
  BAMTYPES="$BAMTYPES,sampleUmiDup,goodCellsUmiDup,annotatedUmiDup,celltypeUmiDup"
  DEDUP="$DEDUP,umi"
fi
OUTDIR="$DATADIR/results_v8"
WORK="$DATADIR/work_v8"
RUNDIR="$DATADIR/nf_run_v8"

if [ "$STAGE" = "run" ]; then
  mkdir -p "$RUNDIR"; cd "$RUNDIR"
  export PS1="${PS1:-}"
  export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  set +u; source ~/.bashrc 2>/dev/null || true; mamba activate scRNAEditing; set -u
  echo "=== $(date) === v8 grid, dataset $DATASET -> $OUTDIR ==="
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
    --dedup "$DEDUP" \
    --strand_modes "annotation,hybrid,caller" \
    --bam_types "$BAMTYPES" \
    --qual 30 \
    --min_cov 10 \
    --bcftools_baq true \
    --requant_editing false \
    --percell_analyses false \
    --mismatch_window_filter true \
    --clustered_filter false \
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
echo "ERROR: unknown STAGE '$STAGE'" >&2; exit 1
