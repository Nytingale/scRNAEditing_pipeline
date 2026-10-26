#!/bin/bash
# Build Replogle K562 CRISPRi-KD ground-truth pseudobulk BAMs for scRNAEditing.
#
# Merges the per-(GEM group, gene) BAMs produced by the a2i_replogle pipeline
# (step 04B, sinto split of MAPQ255 uniquely-mapped reads) across ALL GEM groups
# into one pseudobulk BAM per "arm". This is the KD analog of the teratoma-KO
# merge_gt_bams.sh (which merged T1-T4); here we merge ~200-267 GEM wells.
#
# Source BAMs (CB/UB-tagged, coordinate-sorted, indexed):
#   <SRC>/04B_sample_pert_bams/{KD8_pX_Y}/{GENE}.bam
# Concordance-corrected BAMs (failed-KD-construct cells removed) preferred where
# they exist (only for discordant genes, e.g. ADAR in 29 GEM groups):
#   <SRC>/04B_sample_pert_bams/filtered_bams/{KD8_pX_Y}/{GENE}.bam
#
# Arms (gene -> arm label):
#   ADAR          -> ADAR_KD      (== ADAR1; primary positive control, AEI should drop hard)
#   ADARB1        -> ADARB1_KD    (partial)
#   non-targeting -> NTC_ctrl     (baseline)
# (No ADARB2 arm: ADARB2 is not targeted/expressed in this K562 screen.)
#
# NOTE on reference: Replogle BAMs were aligned with Cell Ranger GRCh38-2020-A
# (UCSC-style chrN naming, GRCh38 primary assembly), coordinate-compatible with
# the scRNAEditing hg38.analysisSet.fa. A quick contig sanity check is printed.
#
# Usage:  qsub scripts/merge_replogle_gt_bams.sh            # all 3 arms
#         qsub scripts/merge_replogle_gt_bams.sh ADAR_KD    # one arm
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -N merge_repl_gt
#$ -pe threaded 4
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
# Do NOT pass -V (SGE JSV / truncated PATH).

set -euo pipefail

export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

THREADS=4

SRC="${EXTERNAL_DIR:-.}/a2i_replogle"
PERTBAMS="${SRC}/04B_sample_pert_bams"
FILTERED="${PERTBAMS}/filtered_bams"

GT="${SCRNA_DIR:-.}/human/replogle"
OUTDIR="${GT}/02_gt_bams"
mkdir -p "$OUTDIR"

# gene:arm pairs
PAIRS=("ADAR:ADAR_KD" "ADARB1:ADARB1_KD" "non-targeting:NTC_ctrl")

# Optional single-arm restriction
WANT_ARM="${1:-}"

for PAIR in "${PAIRS[@]}"; do
    GENE="${PAIR%%:*}"
    ARM="${PAIR##*:}"
    [ -n "$WANT_ARM" ] && [ "$ARM" != "$WANT_ARM" ] && continue

    OUT="${OUTDIR}/merged_${ARM}.bam"
    if [ -f "$OUT" ] && [ -f "${OUT}.bai" ]; then
        echo "=== $(date) === $OUT exists, skipping ==="
        continue
    fi

    # Build file-of-filenames: prefer concordance-corrected BAM, else original.
    FOFN="${OUTDIR}/.fofn_${ARM}_$$"
    : > "$FOFN"
    N_FILT=0; N_ORIG=0
    for D in "${PERTBAMS}"/KD8_*; do          # glob skips the filtered_bams/ subdir
        [ -d "$D" ] || continue
        GEM="$(basename "$D")"
        FSRC="${FILTERED}/${GEM}/${GENE}.bam"
        OSRC="${D}/${GENE}.bam"
        if [ -f "$FSRC" ]; then
            echo "$FSRC" >> "$FOFN"; N_FILT=$((N_FILT+1))
        elif [ -f "$OSRC" ]; then
            echo "$OSRC" >> "$FOFN"; N_ORIG=$((N_ORIG+1))
        fi
    done

    N_IN=$(wc -l < "$FOFN")
    echo "=== $(date) === arm=${ARM} gene='${GENE}'  inputs=${N_IN} (corrected=${N_FILT}, original=${N_ORIG}) ==="
    if [ "$N_IN" -eq 0 ]; then
        echo "ERROR: no source BAMs found for gene '${GENE}'" >&2
        rm -f "$FOFN"
        exit 1
    fi

    samtools merge -@ "$THREADS" -f -b "$FOFN" "$OUT"
    samtools index -@ "$THREADS" "$OUT"
    rm -f "$FOFN"

    NREADS=$(samtools view -@ "$THREADS" -c "$OUT")
    echo "=== $(date) === wrote $OUT  reads=${NREADS} ==="
    echo "    contigs (first 3): $(samtools view -H "$OUT" | awk -F'\t' '/^@SQ/{sub(/SN:/,"",$2);print $2}' | head -3 | tr '\n' ' ')"
    if [ "$NREADS" -eq 0 ]; then
        echo "ERROR: 0 reads in $OUT" >&2
        exit 1
    fi
done

echo "=== $(date) === merge_replogle_gt_bams done ==="
