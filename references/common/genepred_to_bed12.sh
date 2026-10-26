#!/bin/bash
# UCSC genePred table (e.g. ncbiRefSeq.txt.gz, with the leading bin column) -> BED12,
# sorted by chromosome and start. Used for --refseq_bed (RSeQC strandedness inference).
#
# Usage: genepred_to_bed12.sh <ncbiRefSeq.txt.gz> <out.bed>
set -euo pipefail
zcat "$1" | awk 'BEGIN{OFS="\t"} {
    n = $9; split($10, s, ","); split($11, e, ","); sizes = ""; starts = ""
    for (i = 1; i <= n; i++) { sizes = sizes (e[i] - s[i]) ","; starts = starts (s[i] - $5) "," }
    print $3, $5, $6, $2, 0, $4, $7, $8, 0, n, sizes, starts
}' | LC_ALL=C sort -k1,1 -k2,2n > "$2"
echo "$(wc -l < "$2") transcripts -> $2"
