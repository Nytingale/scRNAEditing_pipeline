#!/bin/bash
# Download the built annotation and SNP files from Zenodo into <dir>:
#   <dir>/annotation_files/extended_annots_<genome>.RData   --strand_annots / --alu_annots
#   <dir>/SNP_files/<genome>_SNPs_GRanges.RData              --snp_file
#     (human dbSNP156_common_GRanges.RData, mouse GRCm39_snps_GRanges.RData)
#   <dir>/rediportal/REDIportal_hg38_v2.vcf.gz(.tbi), REDIportal_mm39.vcf   --rediportal
# Usage: bash download_reference_files.sh [dir]    (default: current directory)
set -euo pipefail

VERSION=1.0.0
ZENODO_RECORD=${ZENODO_RECORD:-}          # filled in when the Zenodo record is published
NAME=scRNAEditing_reference_files_${VERSION}
DIR=${1:-.}

[ -n "$ZENODO_RECORD" ] || { echo "Set ZENODO_RECORD (Zenodo record ID; see README.md)" >&2; exit 1; }
URL="https://zenodo.org/records/${ZENODO_RECORD}/files"

mkdir -p "$DIR" && cd "$DIR"
for f in "$NAME.tar.gz" "$NAME.MD5SUMS.txt"; do
    [ -f "$f" ] || curl -fL --retry 3 -o "$f" "$URL/$f?download=1"
done
grep " $NAME.tar.gz\$" "$NAME.MD5SUMS.txt" | md5sum -c -
tar -xzf "$NAME.tar.gz" --strip-components=1
sed "/ $NAME.tar.gz\$/d" "$NAME.MD5SUMS.txt" | md5sum -c --quiet - && echo "All files unpacked and verified."
