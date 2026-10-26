#!/bin/bash
# Symlink SRA FASTQs into CellRanger sample/lane names for T1–T4.
# Pairing: GEO GSE248941 GSM7923327–330 (corrected 2026-09-05).
# CellRanger expects: {Sample}_S{N}_L00{lane}_R{1,2}_001.fastq.gz

set -euo pipefail

FQDIR="${SCRNA_DIR:-.}/human/pm/fastq_files"
CRDIR="${SCRNA_DIR:-.}/human/pm/cellranger_fqs"
mkdir -p "$CRDIR"

# sample  S-index  L001 SRR (large)   L002 SRR (small)
# T1      S1       SRR26992458        SRR26992459
# T2      S2       SRR26992455        SRR26992456
# T3      S3       SRR26992452        SRR26992453
# T4      S4       SRR26992449        SRR26992450

link_lane() {
    local sample="$1" sidx="$2" lane="$3" srr="$4"
    local r
    for r in 1 2; do
        local src="${FQDIR}/${srr}_${r}.fastq.gz"
        local dst="${CRDIR}/${sample}_S${sidx}_L00${lane}_R${r}_001.fastq.gz"
        if [ ! -f "$src" ]; then
            echo "ERROR: missing $src" >&2
            exit 1
        fi
        ln -sfn "$src" "$dst"
        echo "  $dst -> ${srr}_${r}.fastq.gz"
    done
}

echo "=== Remaking CellRanger FASTQ symlinks ==="
link_lane T1 1 1 SRR26992458
link_lane T1 1 2 SRR26992459
link_lane T2 2 1 SRR26992455
link_lane T2 2 2 SRR26992456
link_lane T3 3 1 SRR26992452
link_lane T3 3 2 SRR26992453
link_lane T4 4 1 SRR26992449
link_lane T4 4 2 SRR26992450
echo "=== Done ==="
ls -l "$CRDIR"
