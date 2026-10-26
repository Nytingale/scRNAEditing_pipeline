#!/bin/bash
# Split into per-cell-type BAMs for the celltype-as-sample GT discovery.
#
# METADATA-DRIVEN. This script contains NO dataset constants — no arm list, no
# sample list, no cell-type minimum, no teratoma paths. Everything comes from the
# split plan written by:
#
#   Rscript nextflow/bin/validate_samplesheet.R <sheet.csv> \
#       --dataset <name> --min-cells N --tsv-out publication/CELLTYPE_SCOPE.tsv \
#       --emit-plan <PLANDIR>
#
# so the BAMs emitted here are BY CONSTRUCTION the ones CELLTYPE_SCOPE.tsv
# records as kept. The plan is written only if the sheet validates, and the
# earlier bash re-implementation of the pooled cell counting is gone, so the
# split can no longer drift from the recorded scope.
#
# STRATEGY (matches the original analysis; HANDOFF_2026-09-15 §6):
#   split PER SAMPLE by cell type  ->  merge same cell type across samples
# NOT the reverse. Splitting a merged BAM by barcode misassigns reads, because
# the merged BAM pools samples with raw CB tags, carries no @RG, and 10x barcodes
# are reused across samples. Measured on teratoma: 10 barcodes collide within an
# arm, 8 of them with conflicting cell-type labels. Splitting per sample removes
# that entirely -- each sample's reads are in their own file at the moment the
# label is applied.
#
# BARCODE PREFIXING: after the per-celltype fragments are merged across samples
# the raw CB/CR tags would collide again, and the pipeline dedups AFTER it
# receives the BAM. Dedup is cell-aware (umi_tools --per-cell, MarkDuplicates
# --BARCODE_TAG), so two cells colliding on a barcode look like ONE cell and
# their independent molecules get collapsed. CB and CR are therefore rewritten to
# <dna_code>-<barcode>, which makes keys globally unique and makes
# dedup-after-merge produce exactly what per-sample dedup would. Verified twice:
# 2,587,345 reads round-trip unchanged, 195 @SQ preserved, quickcheck clean, and
# reversing the prefix reproduces the original byte-for-byte.
#
# The code is DNA, not the readable sample_id, because Picard MarkDuplicates
# --BARCODE_TAG rejects anything outside ^[ATCGNatcgn-]*$ -- see the DNA-SAFE
# SAMPLE CODES block below for the encoding and why it is fixed-width base-4.
# `_barcode_prefix_map.tsv` records code -> sample_id so this stays reversible.
# On CEL-Seq2 this is not a nicety -- all 96 of 96 plate barcodes recur across
# all 8 libraries, so (sample_id, barcode) is the ONLY valid cell identity.
#
# TWO TOPOLOGIES (SAMPLESHEET_SPEC.md), branched on `platform`:
#   barcoded (10xv2/10xv3/celseq2) one BAM = many cells -> split by CB, prefix, merge
#   plate    (smartseq2)           one BAM = ONE cell   -> no split, no prefix,
#                                  group whole BAMs by the cell type of cell_id
#
# Writes, per unit (the samplesheet `group`):
#   <OUTBAM>/<unit>__<celltype>.bam(.bai)
#   <OUTBC>/<unit>__<celltype>.tsv    2-col cellKey<TAB>celltype, prefixed to match
#   <FRAGDIR>/_frag_<unit>.csv        samplesheet fragment, incl. the group column
#   <FRAGDIR>/_ok_<unit>              COMPLETION SENTINEL -- see below
#
# THE SENTINEL IS LOAD-BEARING. -hold_jid releases on job COMPLETION, not
# success. Without it, if one unit OOMs while the others succeed the driver would
# build a valid-looking samplesheet and launch the full GT run with that unit
# MISSING -- and for teratoma the NTC control arm is 17 of 43 cell types and
# 9540 of ~14.5k cells, so the benchmark would silently have no control to
# compare against. _ok_<unit> is removed on entry and touched only after `set -e`
# has let the whole task run to completion; the driver refuses to start unless
# every unit in _expected.txt has one.
#
# INPUTS (environment, from submit_celltype_as_sample.sh):
#   PLANDIR  required  directory written by --emit-plan
#   OUTBAM   required  output BAM directory
#   OUTBC    required  output annotation directory
#   FRAGDIR  optional  samplesheet fragments + sentinels (default: PLANDIR)
# The array size is supplied by the submitter as `-t 1-<n units>`; it is
# deliberately NOT hardcoded here.
#$ -N split_ct
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 8
# -l mem, NOT h_vmem: the enforced cgroup cap is min(h_rss, slots * 8 GB) and
# h_vmem is per-thread, so it never set the ceiling. This is a hand-written qsub,
# so the nextflow.config `-l mem` sweep (CHANGELOG 2026-09-15) does NOT cover it.
#
# 16 GB, MEASURED — was 48 GB on a mistaken analogy to SINTO_SPLIT_PERCELL.
# That process splits by CELL, so it opens one output handle per cell and peak
# RAM scales with -p; `filterbarcodes` here groups by CELL TYPE, so it opens one
# handle per cell type (18-23 for teratoma, not thousands). Peak is therefore set
# by `samtools sort -@ N` (~768 MB/thread) plus a handful of handles, and does
# not scale with cell count. Self-test on ADAR1_KO (4 lanes, 639 MB input, 680
# cells, 18 types, job 1703742): MEM_USED 1.58 GB against a 48 GB cap, 3.6 min.
# 16 GB is ~10x the measurement and still above sort's ceiling at 8 threads.
#$ -l mem=16G
#$ -l h_rt=24:00:00

set -eo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

: "${PLANDIR:?PLANDIR not set (run validate_samplesheet.R --emit-plan first)}"
: "${OUTBAM:?OUTBAM not set}"
: "${OUTBC:?OUTBC not set}"
FRAGDIR="${FRAGDIR:-$PLANDIR}"
THREADS="${NSLOTS:-8}"

EXPECTED="${PLANDIR}/_expected.txt"
[ -s "$EXPECTED" ] || { echo "ERROR: missing or empty ${EXPECTED}" >&2; exit 1; }

UNIT=$(sed -n "${SGE_TASK_ID}p" "$EXPECTED")
[ -n "$UNIT" ] || { echo "ERROR: no unit at line ${SGE_TASK_ID} of ${EXPECTED}" >&2; exit 1; }

PLAN="${PLANDIR}/_plan_${UNIT}.tsv"
CELLS="${PLANDIR}/_cells_${UNIT}.tsv"
KEEP="${PLANDIR}/_keep_${UNIT}.txt"
for f in "$PLAN" "$CELLS" "$KEEP"; do
    [ -s "$f" ] || { echo "ERROR: missing or empty plan file ${f}" >&2; exit 1; }
done

mkdir -p "$OUTBAM" "$OUTBC" "$FRAGDIR"

# Clear the sentinel FIRST, so a task that dies part-way can never leave a stale
# success marker from a previous attempt.
OK="${FRAGDIR}/_ok_${UNIT}"
rm -f "$OK"

N_CT=$(wc -l < "$KEEP")
echo "=== $(date) === unit ${UNIT} (task ${SGE_TASK_ID}): ${N_CT} cell types, ${THREADS} threads"

# One topology and one platform per unit: the unit becomes ONE BAM per cell type,
# and barcoded/plate need different handling, so a mixed unit is unbuildable.
TOPOS=$(awk -F'\t' 'NR>1{print $5}' "$PLAN" | sort -u)
PLATFORMS=$(awk -F'\t' 'NR>1{print $4}' "$PLAN" | sort -u)
[ "$(echo "$TOPOS" | wc -l)" -eq 1 ] || {
    echo "ERROR: unit ${UNIT} mixes topologies ($(echo $TOPOS | tr '\n' ' ')); a group cannot span both" >&2; exit 1; }
[ "$(echo "$PLATFORMS" | wc -l)" -eq 1 ] || {
    echo "ERROR: unit ${UNIT} mixes platforms ($(echo $PLATFORMS | tr '\n' ' '))" >&2; exit 1; }
TOPOLOGY="$TOPOS"
PLATFORM="$PLATFORMS"
echo "  platform=${PLATFORM} topology=${TOPOLOGY}"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# ── DNA-SAFE SAMPLE CODES ────────────────────────────────────────────────────
# >>> THE PREFIX MUST BE VALID DNA. <<<
# Dedup is cell-aware, and Picard MarkDuplicates --BARCODE_TAG validates the tag
# against ^[ATCGNatcgn-]*$ (UmiUtil.getTopStrandNormalizedUmi). A readable
# prefix like `T4_AAVS1_ctrl_` therefore aborts the run with
#   picard.PicardException: UMI found with illegal characters
# which is exactly what killed the 2026-09-17 launch at
# GATK_MARKDUPLICATES (AAVS1_ctrl__Melanoblasts). Underscores, digits and most
# letters are all illegal; only A/C/G/T/N and '-' survive.
#
# So each sample_id is encoded as a FIXED-WIDTH base-4 word over {A,C,G,T}
# (A=0 C=1 G=2 T=3), joined to the barcode with '-':
#   T4_AAVS1_ctrl_GGTAATCAGATCACCT  ->  AAT-GGTAATCAGATCACCT
#
# Why fixed-width base-4 rather than a unary scheme (T1=TA, T2=TAA, ...):
#   * it scales — 20 samples need 3 characters, 64 samples still need 3,
#     whereas unary would need 20 and 64. jzlcel already has 8 samples and
#     jzl10x 3, so the width adapts per dataset.
#   * every code is the same length, so the mapping stays unambiguous even if a
#     separator were ever lost, and the barcode is recoverable by fixed offset.
#   * it is dataset-agnostic: nothing here is teratoma-specific.
#
# Codes are assigned over the sorted unique sample_ids of the WHOLE dataset (all
# per-unit plans), not just this unit, so every array task derives the identical
# map and a code means the same thing everywhere. The map is written out, so the
# transform stays reversible without having to guess the scheme.
PREFIX_MAP="${FRAGDIR}/_barcode_prefix_map.tsv"
ALL_SIDS=$(awk -F'\t' 'FNR>1{print $1}' "${PLANDIR}"/_plan_*.tsv 2>/dev/null | sort -u)
N_SIDS=$(printf '%s\n' "$ALL_SIDS" | grep -c .)
[ "$N_SIDS" -ge 1 ] || { echo "ERROR: no sample_ids found in ${PLANDIR}/_plan_*.tsv" >&2; exit 1; }
CODE_W=1; while [ $((4 ** CODE_W)) -lt "$N_SIDS" ]; do CODE_W=$((CODE_W + 1)); done

# index -> fixed-width base-4 DNA word
dna_code () {
    local n=$1 w=$2 d out=""
    for (( d = 0; d < w; d++ )); do
        case $(( n % 4 )) in 0) out="A${out}";; 1) out="C${out}";; 2) out="G${out}";; 3) out="T${out}";; esac
        n=$(( n / 4 ))
    done
    printf '%s' "$out"
}

# Written by whichever task gets there first; identical content either way, so a
# concurrent write is harmless. Built atomically so a reader never sees a
# half-written map.
if [ ! -s "$PREFIX_MAP" ]; then
    TMPMAP=$(mktemp "${PREFIX_MAP}.XXXXXX")
    printf 'code\tsample_id\n' > "$TMPMAP"
    i=0
    while read -r s; do
        [ -n "$s" ] || continue
        printf '%s\t%s\n' "$(dna_code "$i" "$CODE_W")" "$s" >> "$TMPMAP"
        i=$((i + 1))
    done <<< "$ALL_SIDS"
    mv -f "$TMPMAP" "$PREFIX_MAP"
fi
echo "  barcode prefix: ${N_SIDS} sample(s), width ${CODE_W}, map ${PREFIX_MAP}"

sid_code () { awk -F'\t' -v s="$1" 'NR>1 && $2==s {print $1; exit}' "$PREFIX_MAP"; }

# ── Build per-cell-type fragments ────────────────────────────────────────────
if [ "$TOPOLOGY" = "barcoded" ]; then
    # One row per LANE. Lanes of a sample share a sample_id (the dedup unit) and
    # therefore share its cell set: a cell's reads are spread across its lanes,
    # so every lane BAM is searched for every one of that sample's cells.
    # LANE is only used to keep temp filenames unique.
    LANE=0
    while IFS=$'\t' read -r SID BAM ANNOT PF TOPO CID; do
        LANE=$((LANE + 1))
        [ -s "$BAM" ] || { echo "ERROR: missing BAM ${BAM}" >&2; exit 1; }

        # This sample's kept cells. Barcodes stay RAW here: sinto matches them
        # against the tags in the not-yet-rewritten BAM.
        awk -F'\t' -v s="$SID" 'BEGIN{OFS="\t"} NR>1 && $1==s {print $2, $3}' \
            "$CELLS" > "${WORK}/annot_${LANE}.tsv"

        n=$(wc -l < "${WORK}/annot_${LANE}.tsv")
        if [ "$n" -eq 0 ]; then
            echo "  ${SID} lane ${LANE}: no kept cells, skipping"
            continue
        fi
        echo "  ${SID} lane ${LANE}: splitting ${n} cells"

        SDIR="${WORK}/lane_${LANE}"
        mkdir -p "$SDIR"
        sinto filterbarcodes -b "$BAM" -c "${WORK}/annot_${LANE}.tsv" \
            --outdir "$SDIR" -p "$THREADS"

        # Rewrite CB/CR to <dna_code>-<barcode> in the same pass that sorts, so
        # it costs no extra traversal. Header lines are left untouched.
        SCODE=$(sid_code "$SID")
        [ -n "$SCODE" ] || { echo "ERROR: no DNA code for sample_id ${SID} in ${PREFIX_MAP}" >&2; exit 1; }
        while read -r CT; do
            SRC="${SDIR}/${CT}.bam"
            [ -s "$SRC" ] || continue
            samtools view -h "$SRC" \
                | sed -e "/^@/!s/\tCB:Z:/\tCB:Z:${SCODE}-/g" \
                      -e "/^@/!s/\tCR:Z:/\tCR:Z:${SCODE}-/g" \
                | samtools sort -@ "$THREADS" -o "${WORK}/frag_${CT}__${LANE}.bam" -
            rm -f "$SRC"
        done < "$KEEP"
    done < <(tail -n +2 "$PLAN")
else
    # Plate-based: one BAM IS one cell, named by cell_id. Nothing to split and
    # nothing to prefix -- plate cell_ids are already globally unique and there
    # is no CB tag to rewrite. Each whole BAM is routed to the cell type its
    # cell_id maps to.
    LANE=0
    while IFS=$'\t' read -r SID BAM ANNOT PF TOPO CID; do
        LANE=$((LANE + 1))
        [ -s "$BAM" ] || { echo "ERROR: missing BAM ${BAM}" >&2; exit 1; }
        [ -n "$CID" ] || { echo "ERROR: plate row ${SID} has no cell_id" >&2; exit 1; }

        CT=$(awk -F'\t' -v c="$CID" 'NR>1 && $2==c {print $3; exit}' "$CELLS")
        if [ -z "$CT" ]; then
            echo "  ${SID} (${CID}): not in a kept cell type, skipping"
            continue
        fi
        cp "$BAM" "${WORK}/frag_${CT}__${LANE}.bam"
    done < <(tail -n +2 "$PLAN")
fi

# ── Merge each cell type across samples; emit prefixed annot + samplesheet ───
FRAG="${FRAGDIR}/_frag_${UNIT}.csv"
: > "$FRAG"
N_WROTE=0
while read -r CT; do
    mapfile -t PARTS < <(ls -1 "${WORK}/frag_${CT}__"*.bam 2>/dev/null || true)
    [ "${#PARTS[@]}" -gt 0 ] || { echo "  WARN: no BAM for ${CT}" >&2; continue; }

    SID="${UNIT}__${CT}"
    DST="${OUTBAM}/${SID}.bam"
    if [ "${#PARTS[@]}" -eq 1 ]; then
        mv "${PARTS[0]}" "$DST"
    else
        samtools merge -@ "$THREADS" -f "$DST" "${PARTS[@]}"
        rm -f "${PARTS[@]}"
    fi
    samtools index -@ "$THREADS" "$DST"

    # 2-col annotation matching the tags actually in the merged BAM: prefixed for
    # barcoded data, bare cell_id for plate data.
    BC="${OUTBC}/${SID}.tsv"
    if [ "$TOPOLOGY" = "barcoded" ]; then
        # Must use the SAME DNA code as the BAM rewrite above, or the annotation
        # keys stop matching the tags and every cell silently drops out.
        awk -F'\t' -v ct="$CT" 'BEGIN{OFS="\t"}
             NR==FNR { if (FNR>1) code[$2]=$1; next }
             FNR>1 && $3==ct {
                 if (!($1 in code)) { print "ERROR: no code for " $1 > "/dev/stderr"; exit 1 }
                 print code[$1]"-"$2, $3
             }' "$PREFIX_MAP" "$CELLS" > "$BC"
    else
        awk -F'\t' -v ct="$CT" 'BEGIN{OFS="\t"} NR>1 && $3==ct {print $2, $3}' "$CELLS" > "$BC"
    fi

    echo "${SID},${DST},${BC},${BC},${PLATFORM},\"\",\"\",${UNIT}" >> "$FRAG"
    N_WROTE=$((N_WROTE + 1))
    echo "  wrote ${SID}  ($(samtools view -c "$DST") reads, $(wc -l < "$BC") cells, ${#PARTS[@]} part(s))"
done < "$KEEP"

[ "$N_WROTE" -gt 0 ] || { echo "ERROR: unit ${UNIT} produced no celltype BAM" >&2; exit 1; }

# Only now, after every step above has succeeded under `set -e`.
touch "$OK"
echo "=== $(date) === ${UNIT} done (${N_WROTE} of ${N_CT} celltype BAMs), sentinel ${OK}"
