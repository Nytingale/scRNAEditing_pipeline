#!/bin/bash
# Per-SAMPLE (biological-replicate) teratoma GT discovery.
# Sample level = per-teratoma per-arm: T1..T4 x {ADAR1_KO, ADARB1_KO(ADAR2),
# ADARB2_KO(ADAR3), AAVS1_ctrl, NTC_ctrl} = 20 pseudobulk-per-teratoma BAMs.
# This REPLACES the old merged-arm "sample" tier (which pooled T1-T4). The
# cell-type axis is results_gt_celltype_as_sample/ (each cell type pooled across
# T1-T4, split per teratoma before merging) — so celltype is intentionally NOT
# run here.
#
# bam_types: sample, sampleRmDup, sampleUmiDup, recal, recalq20 (NO celltype).
# callers:   bcftools, reditools3, gatk, jacusa2, scallele.
# strand:    annotation for all; hybrid ONLY for bcftools/gatk (others have no
#            read-strand evidence -> hybrid==annotation, dropped as redundant).
# Picks up the bcftools_mpileup_chrom.sh "<*>" fix (drops non-variant placeholder
# rows at the per-chr streaming stage) via a fresh work-dir (no stale cache).
#$ -N nf_pm_ps
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
RUNDIR="${DATADIR}/nf_run_gt_persample"
mkdir -p "$RUNDIR"
cd "$RUNDIR"

export PS1="${PS1:-}"
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

# Extra tiers (goodCells*, annotated*, celltype*) without editing this file:
#   qsub -v BAM_TYPES=sample,sampleRmDup,sampleUmiDup,recal,recalq20,goodCells,goodCellsRmDup,goodCellsUmiDup <this script>
# GT BAMs already hold only annotated cells, so goodCells/annotated read sets
# match sample.
BAM_TYPES="${BAM_TYPES:-sample,sampleRmDup,sampleUmiDup,recal,recalq20}"

echo "=== $(date) === Launching teratoma GT pipeline (PER-SAMPLE T1-T4) ==="

# >>> PER-CELL IS OFF. params.percell_analyses DEFAULTS TO TRUE. <<<
# SINTO_SPLIT_PERCELL + PERCELL_MPILEUP are memory-heavy, sit on NO
# ground-truth critical path, and keep the Nextflow job in state `r` for many
# hours after every benchmark table is written. On the 09-16 merged run per-cell
# ran 14.7 h past the point where all 5 arms were at 64/64, which made the job's
# own completion signal useless and forced the phase-2 scorer to poll data
# readiness instead. These are caller-comparison runs, so it is disabled.

nextflow run "${NFDIR}/main.nf" \
    -profile sge,conda \
    --input "${DATADIR}/samplesheet_pm_gt_persample.csv" \
    --outdir "${DATADIR}/results_gt_persample" \
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
    -work-dir "${DATADIR}/work_gt_persample" \
    -resume

echo "=== $(date) === Per-sample pipeline complete ==="
