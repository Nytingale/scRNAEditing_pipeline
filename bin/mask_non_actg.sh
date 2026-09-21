#!/usr/bin/env bash
# mask_non_actg.sh — Convert non-ACTG bases in a FASTA to N
#
# Some genome assemblies contain IUPAC ambiguity codes (R, B, Y, S, W, K, M, etc.)
# that can cause false positive variant calls with certain tools. This script
# converts all non-ACTG bases (case-insensitive) to N while preserving headers.
#
# Usage:
#   bash mask_non_actg.sh input.fa > output.fa
#   samtools faidx output.fa
#
# The output FASTA preserves line width and headers from the input.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input.fa> [> output.fa]" >&2
    echo "  Converts non-ACTG bases to N in a FASTA file." >&2
    echo "  Remember to re-index the output: samtools faidx output.fa" >&2
    exit 1
fi

INPUT="$1"

if [[ ! -f "$INPUT" ]]; then
    echo "Error: File not found: $INPUT" >&2
    exit 1
fi

awk '{
    if (/^>/) {
        print
    } else {
        gsub(/[^ACTGactgNn\n]/, "N")
        print
    }
}' "$INPUT"
