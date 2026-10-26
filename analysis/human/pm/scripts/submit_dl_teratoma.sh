#!/bin/bash
# Submit download jobs for all teratoma samples (ADARKO, GSE248941)
# Dailamy et al., 2024. Nature Communications.

WKDIR="${SCRNA_DIR:-.}/human/pm"
LOGDIR="${WKDIR}/logs"
SCRIPT="${WKDIR}/scripts/dl_teratoma.sh"

mkdir -p "$LOGDIR"
mkdir -p "${WKDIR}/fastq_files"

while IFS=$'\t' read -r srrid sname teratoma; do
    shortname="${srrid: -3}"
    echo "Submitting ${srrid} (${sname}, ${teratoma})..."
    qsub -N "dl${shortname}" \
         -v WKDIR="${WKDIR}",SRRID="${srrid}" \
         -o "${LOGDIR}/dl_${srrid}.log" \
         "$SCRIPT"
done < "${WKDIR}/hmn_pm_name_match.txt"

echo "All jobs submitted. Check queue with: qstat -u \$USER"
