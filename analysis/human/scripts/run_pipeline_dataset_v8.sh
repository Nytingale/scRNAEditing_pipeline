#!/bin/bash
# =============================================================================
# run_pipeline_dataset_v8.sh  DATASET=<name>  STAGE=run|collect   (2026-10-18)
#
# Full configuration grid on a dataset with the teratoma v8 settings (human/pm/scripts/run_pipeline_pm_gt_v8.sh),
# for ranking configurations by the A-to-G composite rank (every dataset, now and future, goes
# through all configurations). Add a dataset to the case block below.
#   jzlxv3 : Ding et al. PBMC1, 10x 3' v3 (4 libraries)              UMI platform
#   jzlss2 : Ding et al. PBMC1 + PBMC2 Smart-seq2, all 384 + 384 sequenced cells (per-cell BAMs; sample = all cells,
#            goodCells via qc_barcodes = the authors' 253 + 273 QC-passing cells; no UMIs: no umi dedup / *UmiDup)
#   jzlcel : Ding et al. PBMC1 + PBMC2 CEL-seq2 (8 libraries), STARsolo BAMs of 01_starsolobams (no template
#            switching, so no TSO clipping)                   UMI platform
#   obim   : Songco-Casey et al. 2022 O. bimaculoides optic lobe, 2 animals, 10x 3' v3 (26-nt barcode read)  UMI platform
# Datasets other than the ground truths (teratoma, Replogle) get every bam type: the sample-level tiers plus goodCells,
# annotated and celltype (each unprocessed, RmDup and, on UMI platforms, UmiDup and UmiMd = UMI-aware MarkDuplicates,
# added 2026-11-09; -resume reuses every other task). Cell-type call sets are published
# under <sample>__<celltype>/.
# FRESH outdir / work-dir / launch dir: <datadir>/{results_v8, work_v8, nf_run_v8}.
# Do not install into the scRNAEditing env while these run (invalidates -resume).
#   qsub -v DATASET=jzlxv3 human/scripts/run_pipeline_dataset_v8.sh
#   qsub -v DATASET=jzlxv3,STAGE=collect human/scripts/run_pipeline_dataset_v8.sh   (after the run)
#   qsub -v DATASET=jzlxv3,UMIMD_ONLY=1,SHEET_OVERRIDE=<sheet with umimd_bam> ...  (-> results_v8_umimd, 2026-11-09)
#   qsub -v DATASET=jzlxv3,UMI_ONLY=1 ...   (seeded umi_tools only -> results_v8_umi; then STAGE=collect with UMI_ONLY=1)
# References and annotations come from a genome params file (GENOME_PARAMS, a Nextflow -params-file; one per genome in
# annotations/genomes/, keys = main.nf params: genome_fasta, genome_gtf, main_chroms, snp_file, strand_annots,
# alu_annots, rediportal, known_sites_vcf, refseq_bed, star_index; omit what a genome lacks). A dataset not in the case
# block runs with DATADIR, SHEET, UMI and GENOME_PARAMS given directly:
#   qsub -v DATASET=new,DATADIR=<dir>,SHEET=<csv>,UMI=1,GENOME_PARAMS=annotations/genomes/<genome>.yaml ...
# =============================================================================
#$ -N nf_ds_v8
# -q <queue>: add as a "#$" line if your cluster needs one
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
GENOMES="${BASEDIR}/annotations/genomes"
DATASET=${DATASET:?DATASET=jzlxv3|jzlss2|jzlcel|obim or DATADIR + SHEET + UMI + GENOME_PARAMS}
STAGE=${STAGE:-run}

# Defaults per dataset; DATADIR / SHEET / UMI / GENOME_PARAMS from the environment take precedence
case "$DATASET" in
  jzlxv3) d="$BASEDIR/human/jzlxv3";   sh="$d/samplesheet_jzlxv3.csv";         u=1; g="$GENOMES/hg38_analysisSet.yaml" ;;
  jzlss2) d="$BASEDIR/human/jzlss2";   sh="$d/samplesheet_jzlss2_all768.csv";  u=0; g="$GENOMES/hg38_analysisSet.yaml" ;;
  jzlcel) d="$BASEDIR/human/jzlcel";   sh="$d/samplesheet_jzlcel.csv";         u=1; g="$GENOMES/hg38_analysisSet.yaml" ;;
  obim)   d="$BASEDIR/octopus/cmn";    sh="$d/samplesheet_obim.csv";           u=1; g="$GENOMES/obim_ASM119413v2.yaml" ;;
  *)      d=""; sh=""; u=""; g="" ;;
esac
DATADIR=${DATADIR:-$d}; SHEET=${SHEET:-$sh}; UMI=${UMI:-$u}; GENOME_PARAMS=${GENOME_PARAMS:-$g}
for v in DATADIR SHEET UMI GENOME_PARAMS; do
  [ -n "${!v}" ] || { echo "ERROR: DATASET '$DATASET' is not in the case block; set $v" >&2; exit 1; }
done
[ -s "$GENOME_PARAMS" ] || { echo "ERROR: genome params file '$GENOME_PARAMS' not found" >&2; exit 1; }
GENOME_PARAMS=$(readlink -f "$GENOME_PARAMS")
BAMTYPES="sample,sampleRmDup,recal,recalq20,goodCells,goodCellsRmDup,annotated,annotatedRmDup,celltype,celltypeRmDup"
# BQSR needs known sites
grep -q '^known_sites_vcf:' "$GENOME_PARAMS" || BAMTYPES=$(echo "$BAMTYPES" | sed 's/,recal,recalq20//')
DEDUP="none,coordinate"
if [ "$UMI" -eq 1 ]; then
  BAMTYPES="$BAMTYPES,sampleUmiDup,goodCellsUmiDup,annotatedUmiDup,celltypeUmiDup"
  BAMTYPES="$BAMTYPES,sampleUmiMd,goodCellsUmiMd,annotatedUmiMd,celltypeUmiMd"
  DEDUP="$DEDUP,umi,umimd"
fi
SHEET=${SHEET_OVERRIDE:-$SHEET}
OUTDIR="$DATADIR/results_v8"
WORK="$DATADIR/work_v8"
RUNDIR="$DATADIR/nf_run_v8"
# UMIMD_ONLY=1: only the four *UmiMd bam types, into separate dirs; for datasets whose results_v8 predates the current
# CELLSNP_LITE / INFER_STRANDEDNESS (cache miss would rerun and republish every bam type)
if [ "${UMIMD_ONLY:-0}" = 1 ]; then
  [ "$UMI" -eq 1 ] || { echo "ERROR: $DATASET has no UMIs" >&2; exit 1; }
  BAMTYPES="sampleUmiMd,goodCellsUmiMd,annotatedUmiMd,celltypeUmiMd"; DEDUP="umimd"
  OUTDIR="${OUTDIR}_umimd"; WORK="${WORK}_umimd"; RUNDIR="${RUNDIR}_umimd"
fi
# UMI_ONLY=1 (v10): only the four *UmiDup bam types (umi_tools --random-seed=params.umi_seed), into
# separate dirs, replacing the unseeded umi_tools outputs of results_v8 (benchmark only; other analyses use *UmiMd)
if [ "${UMI_ONLY:-0}" = 1 ]; then
  [ "$UMI" -eq 1 ] && [ "${UMIMD_ONLY:-0}" != 1 ] || { echo "ERROR: UMI_ONLY needs a UMI dataset, not with UMIMD_ONLY" >&2; exit 1; }
  BAMTYPES="sampleUmiDup,goodCellsUmiDup,annotatedUmiDup,celltypeUmiDup"; DEDUP="umi"
  OUTDIR="${OUTDIR}_umi"; WORK="${WORK}_umi"; RUNDIR="${RUNDIR}_umi"
fi

if [ "$STAGE" = "run" ]; then
  mkdir -p "$RUNDIR"; cd "$RUNDIR"
  export PS1="${PS1:-}"
  export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  set +u; source ~/.bashrc 2>/dev/null || true; mamba activate scRNAEditing; set -u
  # Slurm since 2026-11-14 (SGE cluster retired)
  PROFILE=(-profile sge,conda)
  [ -n "${SLURM_JOB_ID:-}" ] && PROFILE=(-profile slurm,conda -c "${NFDIR}/conf/slurm.config")
  echo "=== $(date) === v8 grid, dataset $DATASET ($(basename "$GENOME_PARAMS")) -> $OUTDIR ==="
  nextflow run "${NFDIR}/main.nf" \
    "${PROFILE[@]}" \
    -params-file "$GENOME_PARAMS" \
    --input "$SHEET" \
    --outdir "$OUTDIR" \
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
  IDX=$(mktemp -p "$BASEDIR/tmp"); trap 'rm -f "$IDX"' EXIT
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
