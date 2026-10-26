#!/bin/bash
# REDIportal v2 (hg38) TABLE1 -> minimal VCF for --rediportal (overlap summary) and
# STAR --varVCFfile (WASP). A-to-I is A>G on the plus strand and T>C on the minus strand.
#
# Usage: REFDIR=<dir> bash 04_rediportal_vcf.sh
# Output: $REFDIR/REDIportal_hg38_v2.vcf(.gz, .gz.tbi)
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
mkdir -p "$REFDIR" && cd "$REFDIR"

URL="http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v2.txt.gz"
[ -f TABLE1_hg38_v2.txt.gz ] || wget -c -O TABLE1_hg38_v2.txt.gz "$URL"

# TABLE1 columns: 1 Region, 2 Position (1-based), 3 Ref, 4 Ed, 5 Strand, ...
{
    echo '##fileformat=VCFv4.2'
    echo '##source=REDIportal_v2_hg38'
    echo '##INFO=<ID=ED,Number=1,Type=String,Description="RNA editing type">'
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    zcat TABLE1_hg38_v2.txt.gz | awk -F'\t' 'NR > 1 {
        ref = toupper($3)
        if (ref == "A") alt = "G"; else if (ref == "T") alt = "C"; else next
        printf "%s\t%s\t.\t%s\t%s\t.\tPASS\tED=AI\tGT\t0/1\n", $1, $2, ref, alt
    }' | sort -k1,1V -k2,2n
} > REDIportal_hg38_v2.vcf

bgzip -c REDIportal_hg38_v2.vcf > REDIportal_hg38_v2.vcf.gz
tabix -f -p vcf REDIportal_hg38_v2.vcf.gz
echo "$(grep -vc '^#' REDIportal_hg38_v2.vcf) sites -> $REFDIR/REDIportal_hg38_v2.vcf.gz"
