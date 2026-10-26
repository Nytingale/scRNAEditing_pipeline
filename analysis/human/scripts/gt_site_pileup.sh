#!/bin/bash
# gt_site_pileup.sh  (dataset-agnostic GT evaluator, step 2/3 — preferred)
#
# Independent re-quantification of editing at discovered sites using the
# project's established stranded pileup parser (misc_scripts/parse_pileup_stranded.pl).
# This is the strand-aware analog of gt_site_pileup.py and keeps provenance with
# the rest of the codebase.
#
#   samtools mpileup -B -q MQ -Q BQ -d DEPTH -l sites.bed -f ref.fa BAM
#     | parse_pileup_stranded.pl <forward|reverse|unstranded>
#
# 10x 3' libraries are forward (second-strand) -> libtype "forward".
# -B disables BAQ (recommended for RNA editing). Positions with zero coverage
# are simply absent from the output (treated as cov 0 downstream).
#
# Output columns (from the perl): chr pos ref coverage variant_cov A C G T a c g t
#                                 ref_count_plus ref_count_minus editing_level
#
# Usage:
#   gt_site_pileup.sh --bam BAM --bed union_sites.bed --fasta ref.fa --out pileup_ARM.tsv \
#       [--libtype forward] [--min-bq 30] [--min-mq 30] [--max-depth 1000000]

set -euo pipefail

BAM="" BED="" FASTA="" OUT="" LIBTYPE="forward" MINBQ=30 MINMQ=30 MAXDEPTH=1000000
PERL_SCRIPT="${SCRNA_DIR:-.}/misc_scripts/parse_pileup_stranded.pl"

while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAM="$2"; shift 2;;
    --bed) BED="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --libtype) LIBTYPE="$2"; shift 2;;
    --min-bq) MINBQ="$2"; shift 2;;
    --min-mq) MINMQ="$2"; shift 2;;
    --max-depth) MAXDEPTH="$2"; shift 2;;
    --perl) PERL_SCRIPT="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

for v in BAM BED FASTA OUT; do
  [ -n "${!v}" ] || { echo "ERROR: --${v,,} required" >&2; exit 1; }
done
[ -s "$BAM" ]   || { echo "ERROR: missing BAM $BAM" >&2; exit 1; }
[ -s "$BED" ]   || { echo "ERROR: missing BED $BED" >&2; exit 1; }
[ -s "$FASTA" ] || { echo "ERROR: missing FASTA $FASTA" >&2; exit 1; }
[ -f "$PERL_SCRIPT" ] || { echo "ERROR: missing perl parser $PERL_SCRIPT" >&2; exit 1; }
[ -f "${FASTA}.fai" ] || samtools faidx "$FASTA"

echo "=== $(date) === mpileup+parse: $(basename "$BAM")  libtype=${LIBTYPE}  BQ>=${MINBQ} MQ>=${MINMQ} ==="
samtools mpileup -B -q "$MINMQ" -Q "$MINBQ" -d "$MAXDEPTH" \
    -l "$BED" -f "$FASTA" "$BAM" 2>/dev/null \
  | perl "$PERL_SCRIPT" "$LIBTYPE" > "$OUT"

echo "=== $(date) === wrote ${OUT}  ($(($(wc -l < "$OUT") - 1)) covered positions) ==="
