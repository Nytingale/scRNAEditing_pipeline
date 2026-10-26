#!/bin/bash
# Download human + mouse genomes and GTFs for CellRanger mkref
# Source: Ensembl release 98 (matches CellRanger recommended)
# Date: 2026-09-02

set -euo pipefail

REFDIR="${SCRNA_DIR:-.}/human/pm/cellranger_ref"
mkdir -p "$REFDIR"
cd "$REFDIR"

echo "=== Downloading mouse genome (GRCm38) from Ensembl r98 ==="
wget -c https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz

echo "=== Downloading human genome (GRCh38) from Ensembl r98 ==="
wget -c https://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz

echo "=== Decompressing (keeping .gz copies) ==="
for f in Mus_musculus.GRCm38.dna.primary_assembly.fa.gz Mus_musculus.GRCm38.98.gtf.gz \
         Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.98.gtf.gz; do
  zcat "$f" > "${f%.gz}" && echo "  decompressed: $f"
done

echo "=== Recording provenance ==="
cat > "$REFDIR/download_provenance.txt" << EOF
Downloaded: $(date '+%Y-%m-%d %H:%M:%S')
Source: Ensembl release 98
Human: Homo_sapiens.GRCh38.dna.primary_assembly.fa + GRCh38.98.gtf
Mouse: Mus_musculus.GRCm38.dna.primary_assembly.fa + GRCm38.98.gtf
Purpose: CellRanger mkref dual-species reference for teratoma analysis
EOF
ls -lh "$REFDIR/" >> "$REFDIR/download_provenance.txt"

echo "=== Done ==="
