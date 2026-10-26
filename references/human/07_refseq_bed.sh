#!/bin/bash
# hg38_RefSeq.bed (--refseq_bed, RSeQC infer_experiment strandedness): the BED12 distributed
# with RSeQC (file dated 2021-12-13, 172,809 transcripts, md5 8c83e9591e1156ffa56a05bb94dd12ae
# after gunzip).
#
# Usage: REFDIR=<dir> bash 07_refseq_bed.sh
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR" && cd "$REFDIR"
wget -O hg38_RefSeq.bed.gz 'https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz/download'
gunzip -f hg38_RefSeq.bed.gz
echo "md5 $(md5sum < hg38_RefSeq.bed | cut -d' ' -f1) (expected 8c83e9591e1156ffa56a05bb94dd12ae)"
