#!/bin/bash
# mm39_RefSeq.bed (--refseq_bed): BED12 from the UCSC mm39 ncbiRefSeq table.
#
# Usage: REFDIR=<dir> bash 06_refseq_bed.sh
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$REFDIR" && cd "$REFDIR"
[ -f ncbiRefSeq.txt.gz ] || wget -O ncbiRefSeq.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/ncbiRefSeq.txt.gz
bash "$HERE/../common/genepred_to_bed12.sh" ncbiRefSeq.txt.gz mm39_RefSeq.bed
