#!/bin/bash
# =============================================================================
# run_pipeline_aorta_v8.sh  STAGE=run|collect   (2026-10-26)
#
# Aorta (GSE255661) caller runs on the controls, scored against the aorta truth (truth/eval/), to compare
# configuration ranks with the teratoma. Settings = human/pm/scripts/run_pipeline_pm_gt_v8.sh (teratoma v8)
# with the mm39 resources of mouse/scripts/run_pipeline_mjzl.sh, without --alu_annots (human only; the mm39
# RData has no alu.plus.annot / alu.minus.annot, so DOWNSTREAM_R stopped). One run, 4 samples
# (pipeline_bams/samplesheet_aorta_ctrl.csv): CTRL = merged arm of the 3 control mice (arm tier) and the
# 3 mice (sample tier). Inputs from aorta_pipeline_bams.sh: Myh11+ cells, no dedup, @RG SM = mouse (the
# pipeline's coordinate dedup runs per mouse); umidup_bam = the truth UMI BAMs (as the teratoma).
#   qsub mouse/aorta_gse255661/scripts/run_pipeline_aorta_v8.sh
#   qsub -v STAGE=collect mouse/aorta_gse255661/scripts/run_pipeline_aorta_v8.sh   (after the run)
# =============================================================================
#$ -N nf_aorta_v8
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=64G
#$ -l h_rt=336:00:00
set -euo pipefail
export NXF_OPTS="-Xms1g -Xmx10g"
BASEDIR="${SCRNA_DIR:-.}"
NFDIR="${BASEDIR}/nextflow"
REF="${BASEDIR}/mouse/ref_annots"
DATADIR="${BASEDIR}/mouse/aorta_gse255661"
SHEET="$DATADIR/pipeline_bams/samplesheet_aorta_ctrl.csv"
OUTDIR="$DATADIR/results_v8"
WORK="$DATADIR/work_v8"
RUNDIR="$DATADIR/nf_run_v8"
STAGE=${STAGE:-run}
# UMIMD=1: UMI-aware MarkDuplicates (sampleUmiMd) only, from the prebuilt BAMs of
# human/pm/scripts/umimd_bams.sh (umimd_bam column, samplesheet_aorta_ctrl_umimd.csv), into its own outdir / work dir
# (a -resume of work_v8 is not cached: CELLSNP_LITE and INFER_STRANDEDNESS changed since). Scoring splices.
DEDUP="none,coordinate,umi"; BAMTYPES="sample,sampleRmDup,sampleUmiDup,recal,recalq20"
if [ "${UMIMD:-0}" = 1 ]; then
  SHEET="${SHEET%.csv}_umimd.csv"; DEDUP="umimd"; BAMTYPES="sampleUmiMd"
  OUTDIR="${OUTDIR}_umimd"; WORK="${WORK}_umimd"; RUNDIR="${RUNDIR}_umimd"
fi
# UMI=1 (v10): umi_tools configurations only (sampleUmiDup) on the seed-8787 truth/bams, into
# results_v8_umi; STAGE=collect then runs STAGE=place (human/pm/scripts/place_umi_calls.sh: sampleUmiDup subtrees into
# results_v8, Cellsnp-lite hybrid into gt_eval_v8/cellsnp_hybrid/results/results_v8/ as cellsnp_hybrid_v8.sh).
#   qsub -v UMI=1 mouse/aorta_gse255661/scripts/run_pipeline_aorta_v8.sh
PUBDIR=$OUTDIR
if [ "${UMI:-0}" = 1 ]; then
  DEDUP="umi"; BAMTYPES="sampleUmiDup"
  OUTDIR="${OUTDIR}_umi"; WORK="${WORK}_umi"; RUNDIR="${RUNDIR}_umi"
fi

if [ "$STAGE" = "run" ]; then
  [ -s "$SHEET" ] || { echo "ERROR: missing $SHEET (aorta_pipeline_bams.sh merge)" >&2; exit 1; }
  mkdir -p "$RUNDIR"; cd "$RUNDIR"
  export PS1="${PS1:-}"
  export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  set +u; source ~/.bashrc 2>/dev/null || true; mamba activate scRNAEditing; set -u
  echo "=== $(date) === aorta v8 caller run -> $OUTDIR ==="
  nextflow run "${NFDIR}/main.nf" \
    -profile sge,conda \
    --input "$SHEET" \
    --outdir "$OUTDIR" \
    --genome_fasta "${REF}/mm39.fa" \
    --genome_gtf "${REF}/mm39.ncbiRefSeq.gtf" \
    --snp_file "${REF}/GRCm39_snps_GRanges.RData" \
    --strand_annots "${BASEDIR}/annotations/annotation_files/extended_annots_mm39.RData" \
    --rediportal "${REF}/REDIportal_mm39.vcf" \
    --known_sites_vcf "${REF}/GRCm39_snps.vcf.gz" \
    --refseq_bed "${REF}/mm39_RefSeq.bed" \
    --main_chroms "${REF}/chrlist_main_chrs.txt" \
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
    --save_intermediates 4 \
    -work-dir "$WORK" \
    -resume
  echo "=== $(date) === Pipeline complete ==="
  [ "${UMI:-0}" = 1 ] && /cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub -v UMI=1,STAGE=collect "$DATADIR/scripts/run_pipeline_aorta_v8.sh"
  exit 0
fi

if [ "$STAGE" = "place" ]; then
  [ "${UMI:-0}" = 1 ] || { echo "ERROR: STAGE=place needs UMI=1" >&2; exit 1; }
  PLACE_PUB="$PUBDIR" PLACE_HYB="$DATADIR/gt_eval_v8/cellsnp_hybrid/results/results_v8" \
    bash "$BASEDIR/human/pm/scripts/place_umi_calls.sh" "$OUTDIR"
  exit $?
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
  [ "${UMI:-0}" = 1 ] && { STAGE=place bash "$DATADIR/scripts/run_pipeline_aorta_v8.sh"; exit $?; }
  exit 0
fi
echo "ERROR: unknown STAGE '$STAGE'" >&2; exit 1
