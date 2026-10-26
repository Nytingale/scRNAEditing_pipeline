#!/bin/bash
# Download the figure inputs and source data from Zenodo and unpack them into
# this analysis/ directory, at the paths the plot scripts read.
#   bash download_source_data.sh
set -euo pipefail

VERSION=1.0.0
ZENODO_RECORD=${ZENODO_RECORD:-}          # filled in when the Zenodo record is published
NAME=scRNAEditing_analysis_data_${VERSION}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[ -n "$ZENODO_RECORD" ] || { echo "Set ZENODO_RECORD (Zenodo record ID; see README.md)" >&2; exit 1; }
URL="https://zenodo.org/records/${ZENODO_RECORD}/files"

cd "$HERE"
for f in "$NAME.tar.gz" "$NAME.MD5SUMS.txt"; do
    [ -f "$f" ] || curl -fL --retry 3 -o "$f" "$URL/$f?download=1"
done
grep " $NAME.tar.gz\$" "$NAME.MD5SUMS.txt" | md5sum -c -
tar -xzf "$NAME.tar.gz" --strip-components=1
sed "/ $NAME.tar.gz\$/d" "$NAME.MD5SUMS.txt" | md5sum -c --quiet - && echo "All files unpacked and verified."
# Some plot scripts also save a copy of their figure here and do not create the folder.
mkdir -p review/response_figures_2026-10-15
