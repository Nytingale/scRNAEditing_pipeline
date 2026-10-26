#!/bin/bash
# ============================================================================
# WASP_dataset_config.sh — resolve a dataset name to its WASP paths/parameters
# ============================================================================
# Sourced by WASP_prep_pileups.sh and WASP_filter_bams.sh. Reads the SAME
# registry that WASP_common.R reads (WASP_datasets.tsv), so a dataset is
# described in exactly one place and bash and R cannot disagree about it.
#
#   source WASP_dataset_config.sh
#   wasp_load_dataset jzlss2
#
# Sets: DATASET LABEL SPECIES PLATFORM TOPOLOGY LIBTYPE ALIGNER MINCOV
#       DATA_DIR ALLDIR PASSDIR ANALYSIS_DIR REF CHRLIST REDI SRR_LIST
#
# LIBTYPE is the one parameter that is easy to get silently wrong: it is the
# pileup strandedness, and it is a property of the library chemistry, not of
# this analysis. 10x/CEL-Seq2 3' libraries are second-strand ("forward");
# Smart-seq2 is UNSTRANDED. Using "forward" on Smart-seq2 would discard
# roughly half the reads at every site and quietly halve the coverage.
# ============================================================================

wasp_load_dataset() {
    local ds="${1:?usage: wasp_load_dataset <dataset>}"
    local here; here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local registry="${here}/WASP_datasets.tsv"
    local repo; repo="$(cd "${here}/../../.." && pwd)"

    [ -s "$registry" ] || { echo "ERROR: missing registry $registry" >&2; return 1; }

    local row
    row=$(awk -F'\t' -v d="$ds" 'NR>1 && $1==d {print; found=1} END{exit !found}' \
              "$registry") || {
        echo "ERROR: unknown dataset '${ds}'. Known:" >&2
        awk -F'\t' 'NR>1 {print "  " $1}' "$registry" >&2
        return 1
    }

    DATASET="$ds"
    LABEL=$(echo "$row"    | cut -f2)
    SPECIES=$(echo "$row"  | cut -f3)
    PLATFORM=$(echo "$row" | cut -f4)
    TOPOLOGY=$(echo "$row" | cut -f5)
    LIBTYPE=$(echo "$row"  | cut -f6)
    ALIGNER=$(echo "$row"  | cut -f7)
    MINCOV=$(echo "$row"   | cut -f8)

    DATA_DIR="${repo}/human/${ds}"
    ALLDIR="${DATA_DIR}/01_starsolobams_wasp"
    PASSDIR="${DATA_DIR}/01_starsolobams_wasp_filtered"
    ANALYSIS_DIR="${DATA_DIR}/wasp_refbias"
    SRR_LIST="${DATA_DIR}/srr_array_list.txt"

    case "$SPECIES" in
        human)
            REF="${repo}/human/ref_annots/hg38.analysisSet.fa"
            CHRLIST="${repo}/human/ref_annots/chrlist_main_chrs.txt"
            REDI="${repo}/human/ref_annots/REDIportal_hg38_v2.vcf.gz"
            ;;
        *)
            # Mouse is deliberately not wired up: mouse/ref_annots has MGP SNPs
            # but NO REDIportal table, and this diagnostic is REDIportal-guided.
            echo "ERROR: species '${SPECIES}' has no REDIportal table in this repo" >&2
            return 1
            ;;
    esac

    export DATASET LABEL SPECIES PLATFORM TOPOLOGY LIBTYPE ALIGNER MINCOV
    export DATA_DIR ALLDIR PASSDIR ANALYSIS_DIR SRR_LIST REF CHRLIST REDI

    cat <<EOF
--- dataset ${DATASET} (${LABEL}) ---
  platform   : ${PLATFORM}   topology: ${TOPOLOGY}   aligner: ${ALIGNER}
  pileup     : libtype=${LIBTYPE}  mincov=${MINCOV}
  WASP BAMs  : ${ALLDIR}
  filtered   : ${PASSDIR}
  analysis   : ${ANALYSIS_DIR}
EOF
}
