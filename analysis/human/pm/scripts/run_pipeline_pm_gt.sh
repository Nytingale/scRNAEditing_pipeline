#!/bin/bash
# Nextflow on merged CRISPR-arm pseudobulk BAMs (paper-faithful GT).
# AAVS1_ctrl and NTC_ctrl kept as separate samples.
# BAM types: sample / sampleRmDup / recal / celltype.
# celltype is included on purpose: some cell types have very few cells and will
# give poor per-celltype editing calls — that is itself the finding (read depth /
# cell number is a physical limit no calling strategy can overcome).
#$ -N nf_pm_gt
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 1
#$ -l h_vmem=64G

set -eo pipefail
export NXF_OPTS="-Xms1g -Xmx10g"

BASEDIR="${SCRNA_DIR:-.}"
NFDIR="${BASEDIR}/nextflow"
REFDIR="${BASEDIR}/human/ref_annots"
DATADIR="${BASEDIR}/human/pm"
RUNDIR="${DATADIR}/nf_run_gt"
mkdir -p "$RUNDIR"
cd "$RUNDIR"

export PS1="${PS1:-}"
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

echo "=== $(date) === Launching teratoma GT pipeline (merged arms) ==="

nextflow run "${NFDIR}/main.nf" \
    -profile sge,conda \
    --input "${DATADIR}/samplesheet_pm_gt_merged.csv" \
    --outdir "${DATADIR}/results_gt" \
    --genome_fasta "${REFDIR}/hg38.analysisSet.fa" \
    --genome_gtf "${REFDIR}/hg38.ncbiRefSeq.gtf" \
    --snp_file "${BASEDIR}/annotations/SNP_files/dbSNP156_common_GRanges.RData" \
    --strand_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --alu_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --rediportal "${REFDIR}/REDIportal_hg38_v2.vcf.gz" \
    --known_sites_vcf "${REFDIR}/dbSNP156_GRCh38.vcf.gz" \
    --refseq_bed "${BASEDIR}/review/ref_annots/hg38_RefSeq.bed" \
    --main_chroms "${REFDIR}/chrlist_main_chrs.txt" \
    --callers "bcftools,reditools3,gatk,jacusa2" \
    --dedup "none,coordinate" \
    --strand_modes "annotation,hybrid" \
    --bam_types "sample,sampleRmDup,recal,celltype" \
    --qual 30 \
    --min_cov 10 \
    --requant_editing false \
    -work-dir "${DATADIR}/work_gt" \
    -resume

# NOTE: in-pipeline REQUANT_PILEUP is OFF here (it re-piles each BAM per combo,
# redundant + slow). Canonical editing-level columns are added afterwards by the
# pooled post-hoc updater: qsub human/pm/scripts/run_requant_teratoma.sh

echo "=== $(date) === Pipeline complete ==="
