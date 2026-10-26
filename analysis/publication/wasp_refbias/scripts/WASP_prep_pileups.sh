#!/bin/bash
# ============================================================================
# WASP_prep_pileups.sh — Stage 1 source-data generation for the WASP
# reference-bias supplement (Reviewer 2, Major 4). ANY DATASET.
# ============================================================================
# Quantifies the magnitude of STAR reference-mapping bias on A-to-I editing
# WITHOUT re-running the full editing pipeline, by reusing the STAR WASP
# alignments produced by the dataset's stage-0 re-alignment:
#   ALL reads (biased estimate)  : 01_starsolobams_wasp/          (vW-tagged)
#   WASP-passed (bias-corrected) : 01_starsolobams_wasp_filtered/ (vW==1|no vW)
#
# STAR was run with --waspOutputMode SAMtag --varVCFfile REDIportal_hg38_v2.vcf.
# WASP flips the allele at each variant, re-maps the read, and flags reads whose
# alignment is allele-dependent (vW>1) — the reference-biased reads. Comparing
# editing levels at REDIportal A-to-I sites between the two read sets measures
# how much reference bias distorts the editing estimate.
#
# This is a DIAGNOSTIC (human only here — it needs REDIportal). The principled
# cross-species correction is scAllele (reference-free local reassembly),
# already one of the benchmarked callers (Reviewer 2, Major 2).
#
# Steps:
#   1) REDIportal A-to-I sites on main chromosomes -> BED
#   2) merge each set's BAMs -> merged_all.bam / merged_pass.bam
#   3) strand-aware pileup at those sites (engine: human/scripts/gt_site_pileup.sh)
#   4) global read-loss headline counts -> readloss.tsv
# Stage 2 (tables + figure input) is WASP_prep_source_data.R.
#
# Usage:  qsub -N WASP_prep_<ds> publication/wasp_refbias/scripts/WASP_prep_pileups.sh <dataset>
#         bash publication/wasp_refbias/scripts/WASP_prep_pileups.sh <dataset>
# Dataset defaults to jzlxv3. Known datasets: see scripts/WASP_datasets.tsv.
#
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 6
#$ -l mem=60G
#$ -l h_rt=24:00:00
# Do NOT pass -V (SGE JSV / truncated PATH).
# `-l mem` not `-l h_vmem`: this is a hand-written qsub, so the config sweep
# that fixed the scheduler clamp (see CHANGELOG) does not cover it.

set -euo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

THREADS="${NSLOTS:-6}"

# SGE COPIES THE JOB SCRIPT into a spool dir before running it, so under qsub
# ${BASH_SOURCE[0]} is /cm/local/.../job_scripts/<jobid> and the sibling
# WASP_dataset_config.sh is not there. Prefer the real location when running
# interactively, fall back to the canonical path under qsub.
WASP_SCRIPTS_DEFAULT="${SCRNA_DIR:-.}/publication/wasp_refbias/scripts"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/WASP_dataset_config.sh" ] || SCRIPT_DIR="${WASP_SCRIPTS_DEFAULT}"
[ -f "${SCRIPT_DIR}/WASP_dataset_config.sh" ] || {
    echo "ERROR: cannot locate WASP_dataset_config.sh (tried ${SCRIPT_DIR})" >&2; exit 1; }

BASE="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ENG="${BASE}/human/scripts"

# shellcheck source=WASP_dataset_config.sh
source "${SCRIPT_DIR}/WASP_dataset_config.sh"
wasp_load_dataset "${1:-jzlxv3}"

OUT="$ANALYSIS_DIR"
mkdir -p "$OUT"

echo "=== $(date) === WASP reference-bias source data (${DATASET}) ==="

for d in "$ALLDIR" "$PASSDIR"; do
    [ -d "$d" ] || { echo "ERROR: missing BAM dir $d" >&2; exit 1; }
done

# >>> COMPLETENESS GUARD. <<<
# The two arms MUST contain the same units. `-hold_jid` releases on job
# COMPLETION, not success, so if stage 0b dies part-way this job would happily
# merge a partial "pass" set against a complete "all" set and report a WASP
# read loss that is really just missing files — a wrong number that looks
# entirely plausible. Fail loudly instead.
N_ALL=$(find "$ALLDIR"  -maxdepth 1 -name '*Aligned.sortedByCoord.out.bam' | wc -l)
N_PASS=$(find "$PASSDIR" -maxdepth 1 -name '*Aligned.sortedByCoord.out.bam' | wc -l)
if [ "$N_ALL" -ne "$N_PASS" ]; then
    echo "ERROR: unit count mismatch — ${N_ALL} in ${ALLDIR} vs ${N_PASS} in ${PASSDIR}." >&2
    echo "       Stage 0b (WASP_filter_bams.sh ${DATASET}) is incomplete; re-run it." >&2
    exit 1
fi
echo "  units: ${N_ALL} in both arms"

# ── Step 1: REDIportal A-to-I sites on main chromosomes -> BED ──────────────
# Shared across datasets of the same species, so it is built once in the first
# dataset's directory and reused; it depends only on REDI + CHRLIST.
BED="${BASE}/human/ref_annots/reditportal_AI_mainchr.bed"
if [ ! -s "$BED" ]; then
    echo "--- $(date) --- building REDIportal A-to-I BED (main chroms) ---"
    zcat "$REDI" \
      | awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2-1, $2}' \
      | awk 'NR==FNR{keep[$1]=1; next} ($1 in keep)' "$CHRLIST" - \
      | sort -k1,1 -k2,2n > "$BED"
fi
echo "  sites: $(wc -l < "$BED")"

# ── Step 2: merge each set's BAMs (coordinate-sorted; cached) ───────────────
# Plate-based datasets have one BAM PER CELL — 768 for jzlss2 — which blows
# past the open-file limit and the ARG_MAX a glob would need, so the file list
# goes to samtools via -b. Barcoded datasets have a handful of library BAMs
# and take the same path; there is no reason to branch on topology here.
for SET in all pass; do
    if [ "$SET" = "all" ]; then SRC="$ALLDIR"; else SRC="$PASSDIR"; fi
    MB="${OUT}/merged_${SET}.bam"
    LIST="${OUT}/bamlist_${SET}.txt"

    find "$SRC" -maxdepth 1 -name '*Aligned.sortedByCoord.out.bam' \
      | sort > "$LIST"
    N=$(wc -l < "$LIST")
    [ "$N" -ge 1 ] || { echo "ERROR: no BAMs in $SRC" >&2; exit 1; }
    echo "--- ${SET}: ${N} input BAMs ---"

    if [ ! -s "$MB" ]; then
        echo "--- $(date) --- merging ${SET} set from ${SRC} ---"
        # -c -p: plate data has one @RG per cell with colliding IDs; combine
        # identical headers rather than erroring or emitting 768 @RG lines.
        samtools merge -f -@ "$THREADS" -b "$LIST" -c -p "$MB"
        samtools index -@ "$THREADS" "$MB"
    fi
done

# ── Step 3: strand-aware editing-level pileup at REDIportal sites ───────────
# LIBTYPE comes from the registry: 10x/CEL-Seq2 3' libraries are second-strand
# ("forward"), Smart-seq2 is "unstranded". -B (no BAQ), MQ>=30, BQ>=30.
for SET in all pass; do
    PU="${OUT}/pileup_${SET}.tsv"
    if [ ! -s "$PU" ]; then
        echo "--- $(date) --- pileup (${SET}, libtype=${LIBTYPE}) ---"
        bash "${ENG}/gt_site_pileup.sh" \
            --bam "${OUT}/merged_${SET}.bam" \
            --bed "$BED" \
            --fasta "$REF" \
            --out "$PU" \
            --libtype "$LIBTYPE" --min-bq 30 --min-mq 30 --max-depth 1000000
    fi
done

# ── Step 4: read-loss headline (global) ────────────────────────────────────
# NOTE: this is a GENOME-WIDE fraction and is diluted by reads that overlap no
# REDIportal variant (WASP never evaluates those). The loss among evaluated
# reads is computed in Stage 2 from the wasp_filter_summary.csv written by
# WASP_filter_bams.sh, and is much larger. Report both.
echo "--- $(date) --- read-loss counts ---"
A=$(samtools view -c -@ "$THREADS" "${OUT}/merged_all.bam")
P=$(samtools view -c -@ "$THREADS" "${OUT}/merged_pass.bam")
{
    echo -e "metric\tvalue"
    echo -e "reads_all\t${A}"
    echo -e "reads_wasp_pass\t${P}"
    echo -e "reads_removed\t$((A - P))"
    echo -e "frac_removed\t$(awk -v a="$A" -v p="$P" 'BEGIN{ if(a>0) printf "%.6f", (a-p)/a; else print "NA" }')"
} > "${OUT}/readloss.tsv"
cat "${OUT}/readloss.tsv"

echo "=== $(date) === Stage 1 DONE (${DATASET}). Pileups in ${OUT} ==="
echo "Next: Rscript publication/wasp_refbias/scripts/WASP_prep_source_data.R ${DATASET}"
