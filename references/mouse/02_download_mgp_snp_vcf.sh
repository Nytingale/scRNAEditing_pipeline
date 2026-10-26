#!/bin/bash
# ============================================================================
# Download Mouse Genome Project SNP VCF for GATK BQSR and convert chr names
# ============================================================================
#
# Downloads the MGP REL-2112-v8 SNP VCF from EBI, then converts bare
# chromosome names (1, 2, ...) to UCSC-style names (chr1, chr2, ...) for
# compatibility with the mm39 reference.
#
# The raw VCF is ~3.5 GB compressed. bcftools rename takes a while.
#
# Source: EBI Mouse Genomes Project
#   URL: https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/
#   Release: REL-2112-v8 (December 2021)
#   Assembly: GRCm39
#
# Usage:
#   REFDIR=<dir> bash 02_download_mgp_snp_vcf.sh
#
# Output:
#   $REFDIR/GRCm39_snps.vcf.gz       (UCSC chr names, tabix indexed)
#   $REFDIR/GRCm39_snps.vcf.gz.tbi
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

echo "=== $(date) === Working directory: $REFDIR ==="

MGP_URL="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz"
RAW_VCF="mgp_REL2021_snps.vcf.gz"
CONVERTED_VCF="GRCm39_snps.vcf.gz"

# ---- 1. Download raw VCF ----

if [ -f "$CONVERTED_VCF" ]; then
    echo "Converted VCF already exists: $CONVERTED_VCF"
    echo "Delete it to re-download and reconvert."
    exit 0
fi

if [ ! -f "$RAW_VCF" ]; then
    echo "=== $(date) === Downloading MGP SNP VCF (~3.5 GB) ==="
    wget --timestamping "$MGP_URL" -O "$RAW_VCF"
fi

# ---- 2. Build chromosome name mapping ----
# MGP VCF uses bare chromosome names: 1, 2, ..., 19, X, Y, MT
# We need UCSC names: chr1, chr2, ..., chr19, chrX, chrY, chrM

CHRMAP="mgp_to_ucsc_chrmap.txt"

echo "=== $(date) === Building chromosome name mapping ==="
for i in $(seq 1 19) X Y; do
    echo -e "${i}\tchr${i}"
done > "$CHRMAP"
echo -e "MT\tchrM" >> "$CHRMAP"

echo "  Mapped $(wc -l < "$CHRMAP") sequences:"
cat "$CHRMAP"

# ---- 3. Convert chromosome names ----

echo ""
echo "=== $(date) === Converting chromosome names ==="
bcftools annotate \
    --rename-chrs "$CHRMAP" \
    --output-type z \
    --output "$CONVERTED_VCF" \
    --threads 4 \
    "$RAW_VCF"

echo "=== $(date) === Indexing converted VCF ==="
tabix -p vcf "$CONVERTED_VCF"

# ---- 4. Extract SNP coordinates for pipeline RData ----
# The pipeline uses coordinate-based SNP matching for RNA editing site filtering

echo "=== $(date) === Extracting SNP coordinates ==="
bcftools query -f '%CHROM\t%POS\n' "$CONVERTED_VCF" | gzip > GRCm39_snp_coords.gz
echo "  Coordinates: $(zcat GRCm39_snp_coords.gz | wc -l) SNPs"

# ---- 5. Verify ----

echo ""
echo "=== Verification ==="
echo "Converted VCF size: $(ls -lh "$CONVERTED_VCF" | awk '{print $5}')"
echo "Chromosomes in converted VCF:"
bcftools view -h "$CONVERTED_VCF" | grep "^##contig" | head -5
echo "Records per chromosome (sampled):"
for chr in chr1 chr19 chrX; do
    n=$(tabix "$CONVERTED_VCF" "$chr" 2>/dev/null | head -1000 | wc -l)
    echo "  $chr: ${n}+ records (sampled first 1000)"
done

# ---- 6. Provenance ----

cat > download_provenance_mgp_snps.txt <<EOF
# Mouse Genome Project SNP VCF provenance (for GATK BQSR)
# Generated: $(date -u +"%Y-%m-%d %H:%M:%S UTC")

Source VCF:
  File: $RAW_VCF
  Source: EBI Mouse Genomes Project, REL-2112-v8
  URL: $MGP_URL
  Assembly: GRCm39
  Original naming: bare chromosome names (1, 2, ..., 19, X, Y, MT)

Chromosome name conversion:
  Method: bcftools annotate --rename-chrs
  Mapping: bare -> UCSC (1->chr1, ..., MT->chrM)

Converted VCF:
  File: $CONVERTED_VCF
  Naming: UCSC (chr1, chr2, ..., chrM)
  md5: $(md5sum "$CONVERTED_VCF" | awk '{print $1}')

SNP coordinate extract:
  File: GRCm39_snp_coords.gz
  Format: tab-separated CHROM POS
EOF

echo ""
echo "=== $(date) === Done ==="

# Cleanup: remove raw VCF after successful conversion (saves ~3.5 GB)
# Uncomment when satisfied:
# rm -f "$RAW_VCF"
