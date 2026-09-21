#!/bin/bash
# ============================================================================
# Generate scRNAEditing pipeline samplesheet from name match file
# ============================================================================
#
# Name match file format (3-col TSV, no header):
#   SRR_accession    subsample_name    sample_id
#
# Batch match file format (3-col TSV, no header, optional):
#   subsample_name   batch_id    sample_id
#   When provided, batch_id replaces sample_id in the samplesheet and is used
#   for barcode file lookup. This handles cases where a "sample" has multiple
#   GEM wells with reused barcodes (e.g. jzl10x PBMC1 → pbmc1_10x_v2_A/B).
#
# For barcode-based platforms (10xv2, 10xv3, dropseq, celseq2, seqwell, indrop):
#   Looks for barcodes_celltype_{key}.tsv in barcode_dir, where key is:
#     - subsample_name (from name match col 2)
#     - batch_id (from batch match col 2, if provided)
#     - sample_id (from name match col 3)
#
# For Smart-seq2:
#   Uses the celltype merge list (toMerge_celltype_list.txt) as celltype_annot.
#   Each SRR is a separate cell aligned with regular STAR (soloType=None).
#   Only QC-filtered cells (from *_toMerge_qcfiltered_list.txt) become rows.
#
# Usage:
#   bash generate_samplesheet.sh \
#     --name-match <name_match.txt> \
#     --platform <platform> \
#     --bam-dir <bam_directory> \
#     --barcode-dir <barcode_directory> \
#     --output <samplesheet.csv> \
#     [--batch-match <batch_match.txt>]
#
# ============================================================================

set -euo pipefail

NAME_MATCH=""
PLATFORM=""
BAM_DIR=""
BARCODE_DIR=""
OUTPUT=""
BATCH_MATCH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --name-match)   NAME_MATCH="$2"; shift 2;;
        --platform)     PLATFORM="$2"; shift 2;;
        --bam-dir)      BAM_DIR="$2"; shift 2;;
        --barcode-dir)  BARCODE_DIR="$2"; shift 2;;
        --output)       OUTPUT="$2"; shift 2;;
        --batch-match)  BATCH_MATCH="$2"; shift 2;;
        *) echo "Unknown argument: $1"; exit 1;;
    esac
done

if [ -z "$NAME_MATCH" ] || [ -z "$PLATFORM" ] || [ -z "$BAM_DIR" ] || [ -z "$BARCODE_DIR" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: bash generate_samplesheet.sh --name-match <file> --platform <platform> --bam-dir <dir> --barcode-dir <dir> --output <csv> [--batch-match <file>]"
    exit 1
fi

# Load batch match lookup if provided (subsample → batch_id)
declare -A BATCH_MAP
if [ -n "$BATCH_MATCH" ]; then
    while IFS=$'\t' read -r subsample batch_id sample_id_col; do
        BATCH_MAP["$subsample"]="$batch_id"
    done < "$BATCH_MATCH"
    echo "Loaded batch match: ${#BATCH_MAP[@]} entries"
fi

echo "sample_id,bam,celltype_annot,called_barcodes,platform,fastq_r1,fastq_r2,cell_id" > "$OUTPUT"

if [ "$PLATFORM" = "smartseq2" ]; then
    CT_LIST="${BARCODE_DIR}/toMerge_celltype_list.txt"
    if [ ! -f "$CT_LIST" ]; then
        echo "ERROR: celltype merge list not found: $CT_LIST"
        exit 1
    fi

    declare -A QC_PASS
    declare -A QC_FILE
    for qcfile in "${BARCODE_DIR}"/*_toMerge_qcfiltered_list.txt; do
        [ -f "$qcfile" ] || continue
        while IFS=$'\t' read -r subsample sample_id_col; do
            QC_PASS["$subsample"]=1
            QC_FILE["$subsample"]="$qcfile"
        done < "$qcfile"
    done

    while IFS=$'\t' read -r srr subsample sample_id; do
        if [ -n "${QC_PASS[$subsample]+x}" ]; then
            bam="${BAM_DIR}/${srr}_Aligned.sortedByCoord.out.bam"
            qc_list="${QC_FILE[$subsample]}"
            echo "${sample_id},${bam},${CT_LIST},${qc_list},${PLATFORM},,,${subsample}"
        fi
    done < "$NAME_MATCH" >> "$OUTPUT"

else
    while IFS=$'\t' read -r srr subsample sample_id; do
        bam="${BAM_DIR}/${srr}_Aligned.sortedByCoord.out.bam"

        # Determine effective sample_id and barcode lookup key
        effective_id="$sample_id"
        batch_id=""
        if [ -n "${BATCH_MAP[$subsample]+x}" ]; then
            batch_id="${BATCH_MAP[$subsample]}"
            effective_id="$batch_id"
        fi

        # Search for celltype barcode file with fallback chain:
        # batch_id > subsample > sample_id
        ct_file=""
        candidates=("${BARCODE_DIR}/barcodes_celltype_${subsample}.tsv")
        [ -n "$batch_id" ] && candidates=("${BARCODE_DIR}/barcodes_celltype_${batch_id}.tsv" "${candidates[@]}")
        candidates+=("${BARCODE_DIR}/barcodes_celltype_${sample_id}.tsv")

        for candidate in "${candidates[@]}"; do
            if [ -f "$candidate" ]; then
                ct_file="$candidate"
                break
            fi
        done
        if [ -z "$ct_file" ]; then
            echo "WARNING: No celltype barcode file found for ${subsample} (batch=${batch_id}, sample=${sample_id})" >&2
            continue
        fi

        # Search for called_barcodes with same fallback chain
        called_bc=""
        candidates=("${BARCODE_DIR}/barcodes_qcfiltered_${subsample}.tsv")
        [ -n "$batch_id" ] && candidates=("${BARCODE_DIR}/barcodes_qcfiltered_${batch_id}.tsv" "${candidates[@]}")
        candidates+=("${BARCODE_DIR}/barcodes_qcfiltered_${sample_id}.tsv")

        for candidate in "${candidates[@]}"; do
            if [ -f "$candidate" ]; then
                called_bc="$candidate"
                break
            fi
        done

        echo "${effective_id},${bam},${ct_file},${called_bc},${PLATFORM},,,"
    done < "$NAME_MATCH" >> "$OUTPUT"
fi

NROWS=$(tail -n +2 "$OUTPUT" | wc -l)
NSAMPLES=$(tail -n +2 "$OUTPUT" | cut -d',' -f1 | sort -u | wc -l)
echo "Generated: $OUTPUT"
echo "  Rows: $NROWS"
echo "  Samples: $NSAMPLES"
echo "  Platform: $PLATFORM"
