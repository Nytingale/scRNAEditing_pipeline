#!/bin/bash
# Alu elements from the UCSC hg38 RepeatMasker table, split by strand, for the Alu part of
# extended_annots_GRCh38_analysis_set.RData (06_extended_annots.sh). Coordinates stay 0-based
# here; make_extended_annots.R converts them.
#
# Usage: REFDIR=<dir> bash 05_alu_from_rmsk.sh
# Output: $REFDIR/hg38_alu_plus.txt, $REFDIR/hg38_alu_minus.txt
#         (chr, start, end, repName, repFamily, strand)
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR" && cd "$REFDIR"

[ -f rmsk.txt.gz ] || wget -O rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

# rmsk.txt: 6 genoName, 7 genoStart, 8 genoEnd, 10 strand, 11 repName, 13 repFamily
zcat rmsk.txt.gz | grep "Alu" | awk -F"\t" '{if($10 == "+") print $6"\t"$7"\t"$8"\t"$11"\t"$13"\t"$10}' > hg38_alu_plus.txt
zcat rmsk.txt.gz | grep "Alu" | awk -F"\t" '{if($10 == "-") print $6"\t"$7"\t"$8"\t"$11"\t"$13"\t"$10}' > hg38_alu_minus.txt
wc -l hg38_alu_plus.txt hg38_alu_minus.txt
