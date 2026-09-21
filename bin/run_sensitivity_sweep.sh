#!/bin/bash
# ============================================================================
# Post-hoc sensitivity sweep: re-run downstream_filter.R at multiple thresholds
# ============================================================================
# Re-filters existing raw variant calls at different parameter values to show
# that caller rankings and BAM type conclusions are robust to threshold choice.
# Addresses Reviewer 1 Minor 4.
#
# Usage:
#   bash run_sensitivity_sweep.sh <results_dir> <ref_dir> <annot_dir> [ncpus]
#
# Arguments:
#   results_dir  Pipeline results directory (contains <sample>/03_variants/)
#   ref_dir      Reference directory (genome FASTA, GTF)
#   annot_dir    Annotations directory (SNP RData, strand annots RData)
#   ncpus        Number of parallel sweep runs [default: 4]
#
# Output: <results_dir>/sensitivity_sweep/ with collated metrics TSV
# ============================================================================
set -eo pipefail

RESULTS_DIR="$1"
REF_DIR="$2"
ANNOT_DIR="$3"
NCPUS="${4:-4}"

if [ -z "$RESULTS_DIR" ] || [ -z "$REF_DIR" ] || [ -z "$ANNOT_DIR" ]; then
    echo "Usage: bash run_sensitivity_sweep.sh <results_dir> <ref_dir> <annot_dir> [ncpus]"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILTER_R="${SCRIPT_DIR}/downstream_filter.R"
OUTDIR="${RESULTS_DIR}/sensitivity_sweep"
mkdir -p "$OUTDIR"

GENOME_FASTA="${REF_DIR}/hg38.analysisSet.fa"
GENOME_GTF="${REF_DIR}/hg38.ncbiRefSeq.gtf"
SNP_RDATA="${ANNOT_DIR}/SNP_files/dbSNP156_common_GRanges.RData"
STRAND_ANNOTS="${ANNOT_DIR}/annotation_files/extended_annots_GRCh38_analysis_set.RData"
MAIN_CHROMS="${REF_DIR}/chrlist_main_chrs.txt"

# ── Sweep parameter grids (one-at-a-time, others held at default) ────────
declare -A DEFAULTS
DEFAULTS[min_cov]=10
DEFAULTS[min_alt]=3
DEFAULTS[min_edit_level]=0.01
DEFAULTS[jacusa2_min_score]=1.25
DEFAULTS[remove_100pct_editing]=TRUE

declare -A SWEEPS
SWEEPS[min_cov]="5 10 15 20 30 50"
SWEEPS[min_alt]="1 2 3 5 10"
SWEEPS[min_edit_level]="0.005 0.01 0.02 0.05 0.1"
SWEEPS[jacusa2_min_score]="0.5 1.0 1.25 1.5 2.0"
SWEEPS[remove_100pct_editing]="TRUE FALSE"

echo "=== $(date) === Starting sensitivity sweep ==="
echo "Results dir: $RESULTS_DIR"
echo "Output dir:  $OUTDIR"

# ── Find all raw variant files ───────────────────────────────────────────
RAW_FILES=$(find "$RESULTS_DIR" -path "*/03_variants/*" -name "*_variants.txt.gz" -type f 2>/dev/null | sort)
NFILES=$(echo "$RAW_FILES" | wc -l)
echo "Found $NFILES raw variant files"

if [ "$NFILES" -eq 0 ]; then
    echo "ERROR: No raw variant files found in $RESULTS_DIR/*/03_variants/"
    exit 1
fi

# ── Run sweep ────────────────────────────────────────────────────────────
run_filter() {
    local RAW_FILE="$1"
    local PARAM_NAME="$2"
    local PARAM_VALUE="$3"
    local TAG="$4"

    # Parse sample/caller/bam_label from path structure:
    # <results_dir>/<sample>/03_variants/<caller>/<bam_label>/<file>
    local REL_PATH="${RAW_FILE#$RESULTS_DIR/}"
    local SAMPLE=$(echo "$REL_PATH" | cut -d/ -f1)
    local CALLER=$(echo "$REL_PATH" | cut -d/ -f3)
    local BAM_LABEL=$(echo "$REL_PATH" | cut -d/ -f4)

    local SWEEP_OUTDIR="${OUTDIR}/${TAG}/${SAMPLE}/${CALLER}/${BAM_LABEL}"
    mkdir -p "$SWEEP_OUTDIR"
    local METRICS_FILE="${SWEEP_OUTDIR}/metrics.tsv"

    if [ -f "$METRICS_FILE" ]; then
        return 0
    fi

    local MIN_COV=${DEFAULTS[min_cov]}
    local MIN_ALT=${DEFAULTS[min_alt]}
    local MIN_EDIT=${DEFAULTS[min_edit_level]}
    local J2_SCORE=${DEFAULTS[jacusa2_min_score]}
    local RM100=${DEFAULTS[remove_100pct_editing]}

    case "$PARAM_NAME" in
        min_cov)               MIN_COV="$PARAM_VALUE" ;;
        min_alt)               MIN_ALT="$PARAM_VALUE" ;;
        min_edit_level)        MIN_EDIT="$PARAM_VALUE" ;;
        jacusa2_min_score)     J2_SCORE="$PARAM_VALUE" ;;
        remove_100pct_editing) RM100="$PARAM_VALUE" ;;
    esac

    Rscript "$FILTER_R" \
        --input "$RAW_FILE" \
        --caller "$CALLER" \
        --sample "$SAMPLE" \
        --bam_label "$BAM_LABEL" \
        --genome_fasta "$GENOME_FASTA" \
        --snp_rdata "$SNP_RDATA" \
        --strand_annots "$STRAND_ANNOTS" \
        --main_chroms "$MAIN_CHROMS" \
        --min_cov "$MIN_COV" \
        --min_alt "$MIN_ALT" \
        --min_edit_level "$MIN_EDIT" \
        --jacusa2_min_score "$J2_SCORE" \
        --remove_100pct_editing "$RM100" \
        --strand_mode annotation \
        --output_filtered "${SWEEP_OUTDIR}/filtered.txt.gz" \
        --output_summary "${SWEEP_OUTDIR}/summary.txt" \
        --output_stepwise "${SWEEP_OUTDIR}/stepwise.tsv" \
        --output_metrics "$METRICS_FILE" \
        > "${SWEEP_OUTDIR}/log.txt" 2>&1
}

export -f run_filter
export RESULTS_DIR OUTDIR FILTER_R GENOME_FASTA GENOME_GTF SNP_RDATA STRAND_ANNOTS MAIN_CHROMS
export -A DEFAULTS

TOTAL_RUNS=0
for PARAM in min_cov min_alt min_edit_level jacusa2_min_score remove_100pct_editing; do
    for VALUE in ${SWEEPS[$PARAM]}; do
        TAG="${PARAM}_${VALUE}"
        for RAW_FILE in $RAW_FILES; do
            echo "${RAW_FILE}|${PARAM}|${VALUE}|${TAG}"
            ((TOTAL_RUNS++))
        done
    done
done > "${OUTDIR}/_sweep_jobs.txt"

echo "Total sweep runs: $TOTAL_RUNS (using $NCPUS parallel jobs)"

cat "${OUTDIR}/_sweep_jobs.txt" | while IFS='|' read -r RAW_FILE PARAM VALUE TAG; do
    run_filter "$RAW_FILE" "$PARAM" "$VALUE" "$TAG"
done

echo "=== $(date) === Collating sweep metrics ==="

# ── Collate all metrics into one file ────────────────────────────────────
COLLATED="${OUTDIR}/sensitivity_sweep_metrics.tsv"
HEADER_WRITTEN=false

find "$OUTDIR" -name "metrics.tsv" -type f | sort | while read -r MFILE; do
    if [ "$HEADER_WRITTEN" = false ]; then
        # Add sweep_param and sweep_value columns
        SWEEP_TAG=$(echo "$MFILE" | sed "s|${OUTDIR}/||" | cut -d/ -f1)
        SWEEP_PARAM=$(echo "$SWEEP_TAG" | sed 's/_[^_]*$//')
        head -1 "$MFILE" | awk -v OFS='\t' '{print "sweep_param\tsweep_value\t"$0}'
        HEADER_WRITTEN=true
    fi
    SWEEP_TAG=$(echo "$MFILE" | sed "s|${OUTDIR}/||" | cut -d/ -f1)
    # Extract param name and value from tag (e.g., "min_cov_10" -> "min_cov" "10")
    SWEEP_VALUE=$(echo "$SWEEP_TAG" | grep -oP '[^_]+$')
    SWEEP_PARAM=$(echo "$SWEEP_TAG" | sed "s/_${SWEEP_VALUE}$//")
    tail -1 "$MFILE" | awk -v p="$SWEEP_PARAM" -v v="$SWEEP_VALUE" -v OFS='\t' '{print p"\t"v"\t"$0}'
done > "$COLLATED"

NROWS=$(wc -l < "$COLLATED")
echo "Collated $((NROWS - 1)) sweep results to: $COLLATED"
echo "=== $(date) === Sensitivity sweep complete ==="
