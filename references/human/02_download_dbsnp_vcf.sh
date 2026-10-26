#!/bin/bash
# ============================================================================
# Download dbSNP 156 VCF and convert chromosome names for GATK BQSR
# ============================================================================
#
# Downloads dbSNP build 156 VCF for GRCh38 from NCBI, then converts
# RefSeq chromosome accessions (NC_000001.11, etc.) to UCSC-style names
# (chr1, chr2, etc.) for compatibility with our hg38 analysis set reference.
#
# The chromosome name mapping is derived from the NCBI assembly report for
# GCF_000001405.40 (GRCh38.p14), ensuring accuracy.
#
# Source:
#   dbSNP VCF: https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz
#   Assembly report: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
#
# Requirements:
#   - bcftools (for VCF chromosome renaming)
#   - tabix (for VCF indexing)
#   - wget
#
# Usage:
#   REFDIR=<dir> bash 02_download_dbsnp_vcf.sh
#
# Output:
#   $REFDIR/dbSNP156_GRCh38.vcf.gz       (UCSC chr names, indexed)
#   $REFDIR/dbSNP156_GRCh38.vcf.gz.tbi
#
# ============================================================================

set -euo pipefail

REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR"
cd "$REFDIR"

echo "=== $(date) === Working directory: $REFDIR ==="

# ---- 1. Download dbSNP 156 VCF (GRCh38) ----

DBSNP_URL="https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz"
DBSNP_TBI_URL="https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz.tbi"
DBSNP_RAW="GCF_000001405.40.gz"
DBSNP_RAW_TBI="GCF_000001405.40.gz.tbi"
DBSNP_OUT="dbSNP156_GRCh38.vcf.gz"

if [ -f "$DBSNP_OUT" ]; then
    echo "Converted dbSNP VCF already exists: $DBSNP_OUT"
    echo "Delete it to re-download and reconvert."
    exit 0
fi

if [ ! -f "$DBSNP_RAW" ]; then
    echo "=== $(date) === Downloading dbSNP 156 VCF (GRCh38, ~25 GB) ==="
    echo "This may take several hours depending on network speed."
    wget --timestamping "$DBSNP_URL" -O "$DBSNP_RAW"
fi

if [ ! -f "$DBSNP_RAW_TBI" ]; then
    echo "=== $(date) === Downloading tabix index ==="
    wget --timestamping "$DBSNP_TBI_URL" -O "$DBSNP_RAW_TBI"
fi

# ---- 2. Download NCBI assembly report for chromosome name mapping ----

ASSEMBLY_REPORT_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
ASSEMBLY_REPORT="GCF_000001405.40_GRCh38.p14_assembly_report.txt"

if [ ! -f "$ASSEMBLY_REPORT" ]; then
    echo "=== $(date) === Downloading assembly report for chromosome name mapping ==="
    wget "$ASSEMBLY_REPORT_URL" -O "$ASSEMBLY_REPORT"
fi

# ---- 3. Build chromosome name mapping file ----
# Extract: RefSeq accession (col 7) -> UCSC name (col 10)
# Only assembled-molecule and unlocalized/unplaced scaffolds on primary assembly

CHRMAP="ncbi_to_ucsc_chrmap.txt"

echo "=== $(date) === Building chromosome name mapping ==="
awk -F'\t' '
    $0 !~ /^#/ && NF >= 10 && $10 != "na" {
        print $7 "\t" $10
    }
' "$ASSEMBLY_REPORT" > "$CHRMAP"

echo "  Mapped $(wc -l < "$CHRMAP") sequences"
echo "  First 5 entries:"
head -5 "$CHRMAP"

# ---- 4. Convert chromosome names in VCF ----

echo ""
echo "=== $(date) === Converting chromosome names (this takes a while) ==="
bcftools annotate \
    --rename-chrs "$CHRMAP" \
    --output-type z \
    --output "$DBSNP_OUT" \
    --threads 4 \
    "$DBSNP_RAW"

echo "=== $(date) === Indexing converted VCF ==="
tabix -p vcf "$DBSNP_OUT"

# ---- 5. Verify conversion ----

echo ""
echo "=== Verification ==="
echo "Original chromosomes (first 5):"
bcftools view -h "$DBSNP_RAW" | grep "^##contig" | head -5
echo ""
echo "Converted chromosomes (first 5):"
bcftools view -h "$DBSNP_OUT" | grep "^##contig" | head -5
echo ""
echo "Record count (first 3 chrs):"
for chr in chr1 chr2 chr3; do
    n=$(tabix "$DBSNP_OUT" "$chr" 2>/dev/null | head -1000 | wc -l)
    echo "  $chr: ${n}+ records (sampled first 1000)"
done

echo ""
echo "Output: $DBSNP_OUT ($(ls -lh "$DBSNP_OUT" | awk '{print $5}'))"
echo "Index:  ${DBSNP_OUT}.tbi"

# ---- 6. Record provenance ----

cat > download_provenance_dbsnp.txt <<EOF
# dbSNP VCF provenance for GATK BQSR
# Generated: $(date -u +"%Y-%m-%d %H:%M:%S UTC")

Source VCF:
  File: $DBSNP_RAW
  Source: NCBI dbSNP build 156
  URL: $DBSNP_URL
  Assembly: GCF_000001405.40 (GRCh38.p14)
  Original naming: RefSeq accessions (NC_000001.11, etc.)

Chromosome name conversion:
  Mapping source: NCBI assembly report
  Report URL: $ASSEMBLY_REPORT_URL
  Method: bcftools annotate --rename-chrs
  Mapped sequences: $(wc -l < "$CHRMAP")

Converted VCF:
  File: $DBSNP_OUT
  Naming: UCSC (chr1, chr2, etc.)
  md5: $(md5sum "$DBSNP_OUT" | awk '{print $1}')
EOF

echo ""
echo "=== $(date) === Done ==="
echo "Provenance recorded in: download_provenance_dbsnp.txt"

# ---- 7. Optional cleanup ----
# Uncomment to remove the raw NCBI VCF after conversion (saves ~25 GB)
# rm -f "$DBSNP_RAW" "$DBSNP_RAW_TBI" "$ASSEMBLY_REPORT"
