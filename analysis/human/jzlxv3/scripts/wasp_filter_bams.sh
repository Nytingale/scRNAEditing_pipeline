#!/bin/bash
# ============================================================================
# Create WASP-filtered BAMs from WASP-tagged STARsolo output
# ============================================================================
# For each WASP-aligned BAM, creates a filtered version keeping only reads
# that either:
#   (a) do not overlap any variant site (no vW tag), OR
#   (b) passed WASP filtering (vW:i:1)
#
# This removes reference-biased reads (vW:i:2-7) while retaining all reads
# that don't touch editing sites. The result is a bias-corrected BAM.
#
# Also generates a comparison summary of read counts.
#
# Usage: bash wasp_filter_bams.sh
# ============================================================================

set -euo pipefail

WASP_DIR="${SCRNA_DIR:-.}/human/jzlxv3/01_starsolobams_wasp"
FILT_DIR="${SCRNA_DIR:-.}/human/jzlxv3/01_starsolobams_wasp_filtered"
SRRS="SRR9169192 SRR9169193 SRR9169194 SRR9169414"

mkdir -p "$FILT_DIR"

echo "SRR,total_reads,no_variant_overlap,wasp_passed,wasp_failed,pct_removed" > "${FILT_DIR}/wasp_filter_summary.csv"

for SRR in $SRRS; do
    INBAM="${WASP_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"
    OUTBAM="${FILT_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"

    if [ ! -f "$INBAM" ]; then
        echo "WARNING: $INBAM not found, skipping"
        continue
    fi

    echo "=== $(date) === Filtering $SRR ==="

    # Keep reads that don't have vW tag OR have vW:i:1 (WASP-passed)
    # samtools -e expression filter: '![vW] || [vW]==1'
    samtools view -h -b \
        -e '![vW] || [vW]==1' \
        "$INBAM" \
        -o "$OUTBAM"

    samtools index "$OUTBAM"

    # Compute stats
    TOTAL=$(samtools view -c "$INBAM")
    FILTERED=$(samtools view -c "$OUTBAM")
    REMOVED=$((TOTAL - FILTERED))

    WASP_PASS=$(samtools view -c -d vW:1 "$INBAM" 2>/dev/null || echo "0")
    WASP_FAIL=$REMOVED
    NO_VW=$((FILTERED - WASP_PASS))

    if [ "$TOTAL" -gt 0 ]; then
        PCT=$(awk "BEGIN {printf \"%.4f\", $REMOVED/$TOTAL * 100}")
    else
        PCT="0"
    fi

    echo "  Total: $TOTAL | Kept: $FILTERED | Removed: $REMOVED ($PCT%)"
    echo "$SRR,$TOTAL,$NO_VW,$WASP_PASS,$WASP_FAIL,$PCT" >> "${FILT_DIR}/wasp_filter_summary.csv"
done

echo ""
echo "=== $(date) === Filter summary ==="
column -t -s',' "${FILT_DIR}/wasp_filter_summary.csv"

echo ""
echo "Filtered BAMs: $FILT_DIR"
echo "Summary CSV:   ${FILT_DIR}/wasp_filter_summary.csv"
echo ""
echo "Next: create samplesheets for both WASP and WASP-filtered BAMs"
echo "      and run the Nextflow pipeline on each."
