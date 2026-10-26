#!/bin/bash
# =============================================================================
# run_pipeline_pm_gt_v8.sh  TIER=arm|sample|celltype|splitB  STAGE=run|collect   (2026-10-14)
#
# v8 caller runs: all four teratoma caller runs recomputed from scratch, so every caller comes from one
# pipeline version and one conda env (user 2026-10-14). Replaces the v6 runs (results_gt_{v6,persample_v6,
# celltype_v6}, results_gt_split_B), which are kept unchanged except their BCFtools outputs
# (moved to human/pm/_superseded_bcftools_refN_2026-10-14/).
# Why: BCFTOOLS_MPILEUP REF=N race (gt_eval_v8/README.md); a -resume of the v6 work dirs recomputes
# everything anyway (the scRNAEditing conda env changed on 2026-10-10, and the env directory is part of
# every task hash).
# Same inputs and grid as the v6 runs, with the split-B settings for all tiers:
#   --callers + cellsnp           (v6 arm/sample/celltype ran cellSNP-lite separately, cellsnp_v6.sh)
#   --save_intermediates 4        (step-4 / step-5 tables in the work dirs; STAGE=collect copies them
#                                  next to each published pre100 as *_step4_cov_edit_filtered / *_step5_nomm)
#   --mismatch_window_filter true (= downstream_filter.R default, as v6)
#   BCFtools: REF=N race fixed, BAQ ON (--bcftools_baq true, the pipeline default from 2026-10-15)
#   (the 2026-10-14 launch passed --bcftools_baq false, but the string "false" was not coerced, so it ran
#    BAQ on. BAQ on was then adopted (gt_eval_v8/baq_compare/), so those outputs are the final ones.
#    A BAQ-off resume on 2026-10-15 was aborted after ~6 min; see gt_eval_v8/README.md)
# FRESH outdir / work-dir / launch dir per tier: results_gt_<tier>_v8, work_gt_<tier>_v8, nf_run_gt_<tier>_v8
# (arm: results_gt_v8). Do not install into the scRNAEditing env while these run (invalidates -resume).
#   for t in arm sample celltype splitB; do qsub -v TIER=$t human/pm/scripts/run_pipeline_pm_gt_v8.sh; done
#   qsub -v TIER=<tier>,STAGE=collect human/pm/scripts/run_pipeline_pm_gt_v8.sh   (after the run)
# =============================================================================
#$ -N nf_pm_gt_v8
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
DATADIR="${BASEDIR}/human/pm"
TIER=${TIER:?TIER=arm|sample|celltype|splitB}
STAGE=${STAGE:-run}

case "$TIER" in
  arm)      SHEET="$DATADIR/samplesheet_pm_gt_merged_dedupv6.csv";              NAME=gt_v8 ;;
  sample)   SHEET="$DATADIR/samplesheet_pm_gt_persample_dedupv6.csv";           NAME=gt_persample_v8 ;;
  celltype) SHEET="$DATADIR/samplesheet_pm_gt_celltype_as_sample_dedupv6.csv";  NAME=gt_celltype_v8 ;;
  splitB)   SHEET="$DATADIR/gt_eval_v8/split_cells/samplesheet_split_B.csv";    NAME=gt_split_B_v8 ;;
  *) echo "ERROR: unknown TIER '$TIER'" >&2; exit 1 ;;
esac
OUTDIR="$DATADIR/results_$NAME"
WORK="$DATADIR/work_$NAME"
RUNDIR="$DATADIR/nf_run_$NAME"

if [ "$STAGE" = "run" ]; then
  mkdir -p "$RUNDIR"; cd "$RUNDIR"
  export PS1="${PS1:-}"
  export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  set +u; source ~/.bashrc 2>/dev/null || true; mamba activate scRNAEditing; set -u
  echo "=== $(date) === teratoma v8 caller run, tier $TIER -> $OUTDIR ==="
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
    --bcftools_baq true \
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
echo "ERROR: unknown STAGE '$STAGE'" >&2; exit 1
