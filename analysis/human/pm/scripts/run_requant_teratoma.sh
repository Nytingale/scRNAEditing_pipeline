#!/bin/bash
# Post-hoc re-quantification of editing level + coverage for ALL method outputs
# of the teratoma CRISPR-KO GT run. Adds canonical, strand-aware columns
# (edit_level_pileup, cov_pileup, edited_reads, ref_reads) to every
# results_gt/<arm>/04_filtered/<caller>/<bam_type>/<strand>/snp_*/*_{filtered,pre100}.*
# using the MATCHING bam_type BAM (sample->01_merged, sampleRmDup->02_dedup,
# recal->02_gatk_preprocess, ...). Editing level is re-derived directly from the
# reads so values are comparable across callers.
#
# HEAVY (streams multi-GB BAMs) -> MUST be submitted with qsub, never the login node:
#   qsub human/pm/scripts/run_requant_teratoma.sh
#$ -N requant_ter
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 2
#$ -l h_vmem=24G
#$ -l h_rt=24:00:00

set -euo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

BASE="${SCRNA_DIR:-.}"
RESULTS="${BASE}/human/pm/results_gt"
FASTA="${BASE}/human/ref_annots/hg38.analysisSet.fa"

echo "=== $(date) === Requantifying teratoma GT outputs in ${RESULTS} ==="
Rscript "${BASE}/human/scripts/requant_editing.R" \
    --results_dir "${RESULTS}" \
    --fasta       "${FASTA}" \
    --libtype     forward \
    --backup
echo "=== $(date) === Done. Canonical columns added: edit_level_pileup, cov_pileup, edited_reads, ref_reads ==="
