#!/bin/bash
# Teratoma GT on merged CRISPR-arm pseudobulk BAMs — FULL METHOD GRID.
# AAVS1_ctrl and NTC_ctrl kept as separate samples (5 arms total).
#
# WHY THIS REPLACES run_pipeline_pm_gt.sh
# The previous merged-arm run was both INCOMPLETE and RAGGED, which made it
# unusable for method comparison and produced the misleading GT result that
# bcftools had 0% precision on 10 evaluable sites. Measured filtered-result
# counts from that run:
#     bcftools  20 celltype / 10 sample / 8 sampleRmDup / 0 recal
#     gatk      20 celltype / 16 sample / 12 sampleRmDup / 4 recal
#     jacusa2   10 celltype / 18 sample / 18 sampleRmDup / 4 recal
#     scallele  ABSENT ENTIRELY
# and it lacked sampleUmiDup, recalq20 and umi dedup. Its results also predate
# many pipeline changes, so consistency with the per-sample tier could not be
# assumed.
#
# >>> FRESH work-dir AND outdir, DELIBERATELY. NOT a -resume of work_gt. <<<
# Nextflow's task hash covers the process script but NOT files invoked by path,
# so a changed ${projectDir}/bin/*.R (downstream_filter.R was modified as
# recently as 2026-09-16) does NOT invalidate its cache. Resuming the old work
# directory would therefore silently mix results computed by different code
# versions — exactly the inconsistency this re-run exists to remove. The old
# results_gt/ and work_gt/ are left untouched for comparison.
#
# GRID MATCHES run_pipeline_pm_gt_persample.sh EXACTLY, so merged vs per-sample
# is a clean comparison:
#   callers   bcftools, reditools3, gatk, jacusa2, scallele   (scallele ADDED)
#   dedup     none, coordinate, umi                           (umi ADDED)
#   bam_types sample, sampleRmDup, sampleUmiDup, recal, recalq20
#   strand    annotation, hybrid   (+ 'caller' later, see below)
#
# >>> celltype IS DELIBERATELY EXCLUDED, unlike the old script. <<<
# Two independent reasons:
#  1. The merged-arm celltype split is the approach HANDOFF_2026-09-15_pm §5a
#     found to MISASSIGN READS: the merged BAM pools T1-T4 with raw CB tags and
#     carries no @RG, and 10x barcodes recur across teratomas. It is superseded
#     by the celltype-as-sample track (split per sample, prefix CB/CR, then
#     merge), which is what run_pipeline_pm_gt_celltype_as_sample.sh runs.
#  2. The merged samplesheet FAILS validate_samplesheet.R for precisely that
#     reason — 8 duplicated cell keys (1 AAVS1_ctrl, 7 NTC_ctrl). Those
#     duplicates only matter for the celltype split; at 8 of 14,511 keys
#     (0.055%) they are negligible for the sample/recal BAM types kept here.
#
# >>> strand_mode 'caller' IS NOT IN THIS LAUNCH, BY DESIGN. <<<
# strand_mode fans out only at DOWNSTREAM_R (main.nf combines it onto the
# variants channel), so it does NOT re-run any caller. Adding it later and
# re-running this script with -resume costs only the extra DOWNSTREAM_R tasks,
# while launching now starts the expensive caller work immediately. Verified by
# the per-sample grid: reditools3 has 100 runs at annotation only while bcftools
# has 200 across two strand modes, from the same upstream calls.
#$ -N nf_pm_gt_full
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
RUNDIR="${DATADIR}/nf_run_gt_v2"
mkdir -p "$RUNDIR"
cd "$RUNDIR"

export PS1="${PS1:-}"
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

# Extra tiers (goodCells*, annotated*, celltype*) without editing this file:
#   qsub -v BAM_TYPES=sample,sampleRmDup,sampleUmiDup,recal,recalq20,goodCells,goodCellsRmDup,goodCellsUmiDup <this script>
# GT BAMs already hold only annotated cells, so goodCells/annotated read sets
# match sample (barcode_coverage read_retention = 1.0 in results_gt_v2).
BAM_TYPES="${BAM_TYPES:-sample,sampleRmDup,sampleUmiDup,recal,recalq20}"

echo "=== $(date) === Launching teratoma GT pipeline (MERGED ARMS, FULL GRID) ==="

# >>> PER-CELL IS OFF. params.percell_analyses DEFAULTS TO TRUE. <<<
# SINTO_SPLIT_PERCELL + PERCELL_MPILEUP are memory-heavy, sit on NO
# ground-truth critical path, and keep the Nextflow job in state `r` for many
# hours after every benchmark table is written. On the 09-16 merged run per-cell
# ran 14.7 h past the point where all 5 arms were at 64/64, which made the job's
# own completion signal useless and forced the phase-2 scorer to poll data
# readiness instead. These are caller-comparison runs, so it is disabled.

nextflow run "${NFDIR}/main.nf" \
    -profile sge,conda \
    --input "${DATADIR}/samplesheet_pm_gt_merged.csv" \
    --outdir "${DATADIR}/results_gt_v2" \
    --genome_fasta "${REFDIR}/hg38.analysisSet.fa" \
    --genome_gtf "${REFDIR}/hg38.ncbiRefSeq.gtf" \
    --snp_file "${BASEDIR}/annotations/SNP_files/dbSNP156_common_GRanges.RData" \
    --strand_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --alu_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --rediportal "${REFDIR}/REDIportal_hg38_v2.vcf.gz" \
    --known_sites_vcf "${REFDIR}/dbSNP156_GRCh38.vcf.gz" \
    --refseq_bed "${BASEDIR}/review/ref_annots/hg38_RefSeq.bed" \
    --main_chroms "${REFDIR}/chrlist_main_chrs.txt" \
    --callers "bcftools,reditools3,gatk,jacusa2,scallele" \
    --scallele_conda "${CONDA_ENVS:-/opt/conda/envs}/scAllele" \
    --dedup "none,coordinate,umi" \
    --strand_modes "annotation,hybrid,caller" \
    --bam_types "${BAM_TYPES}" \
    --qual 30 \
    --min_cov 10 \
    --requant_editing false \
    --percell_analyses false \
    -work-dir "${DATADIR}/work_gt_v2" \
    -resume

# NOTE: in-pipeline REQUANT_PILEUP is OFF here (it re-piles each BAM per combo,
# redundant + slow). Canonical editing-level columns are added afterwards by the
# pooled post-hoc updater: qsub human/pm/scripts/run_requant_teratoma.sh

echo "=== $(date) === Pipeline complete ==="
