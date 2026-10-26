#!/bin/bash
# Split one teratoma merged BAM to one CRISPR-arm barcode group.
# Usage: bash split_gt_bams.sh T1 ADAR1_KO
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -N split_gt
#$ -pe threaded 2
#$ -l h_vmem=16G
#$ -l h_rt=8:00:00
# Do NOT pass -V (SGE JSV / truncated PATH).

set -euo pipefail

export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

SAMPLE="${1:?sample T1-T4}"
ARM="${2:?arm e.g. ADAR1_KO}"

BASE="${SCRNA_DIR:-.}/human/pm"
BAM="${BASE}/01_starsolobams/${SAMPLE}_merged.bam"
BC="${BASE}/gt_barcodes/${SAMPLE}_${ARM}.tsv"
OUTDIR="${BASE}/02_gt_bams"
OUT="${OUTDIR}/${SAMPLE}_${ARM}.bam"
TMP="${OUTDIR}/tmp_${SAMPLE}_${ARM}_$$"

mkdir -p "$OUTDIR" "$TMP"

if [ ! -f "$BAM" ]; then
    echo "ERROR: missing $BAM" >&2
    exit 1
fi
if [ ! -s "$BC" ]; then
    echo "ERROR: missing/empty barcode file $BC" >&2
    exit 1
fi

if [ -f "$OUT" ] && [ -f "${OUT}.bai" ]; then
    echo "=== $(date) === $OUT already exists, skipping ==="
    rm -rf "$TMP"
    exit 0
fi

echo "=== $(date) === ${SAMPLE} ${ARM}  barcodes=$(wc -l < "$BC") ==="
# Two-column barcode\tlabel for sinto (col2 unused by filter but required)
awk -F'\t' '{gsub(/[[:space:]]+$/, "", $1); if ($1 != "") print $1"\tfiltered"}' "$BC" \
    | sort -u -k1,1 > "${TMP}/sinto.tsv"

sinto filterbarcodes \
    -b "$BAM" \
    -c "${TMP}/sinto.tsv" \
    --outdir "$TMP" \
    -p 2

# sinto names the output after col2 ("filtered.bam") — rename.
if [ -f "${TMP}/filtered.bam" ]; then
    mv "${TMP}/filtered.bam" "$OUT"
else
    samtools view -h -b "$BAM" -D CB:<(cut -f1 "${TMP}/sinto.tsv") > "$OUT"
fi

samtools index -@ 2 "$OUT"
N=$(samtools view -c "$OUT")
echo "=== $(date) === wrote $OUT  reads=$N ==="
if [ "$N" -eq 0 ]; then
    echo "ERROR: 0 reads retained" >&2
    rm -rf "$TMP"
    exit 1
fi
rm -rf "$TMP"
