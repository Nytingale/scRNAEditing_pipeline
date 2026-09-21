#!/bin/bash
# One-chromosome scAllele worker for SCALLELE (chr-scatter, matches the
# bcftools/gatk pattern). Runs scAllele single-threaded (-n 1) on one contig so
# the module can scatter contigs across the allocated cpus with xargs -P instead
# of the old serial per-chromosome loop. Writes per-chr _snvs_<CHR>.txt (all
# SNVs: PASS+LowQual, for editing discovery) and _bl_<CHR>.txt (PASS-only
# genetic-variant SNP blacklist). Per-chr output files avoid shared-file races.
#
# Args: CHR BAM FASTA MINCOV MINALT MINEDIT QUAL STRAND RUNMODE
set -uo pipefail

CHR="$1"; BAM="$2"; FASTA="$3"; MINCOV="$4"; MINALT="$5"; MINEDIT="$6"
QUAL="$7"; STRAND="$8"; RUNMODE="$9"
PREF="_scallele_${CHR}"
rm -f "${PREF}".*

scAllele \
    -b "$BAM" \
    -g "$FASTA" \
    -o "$PREF" \
    -c "$CHR" \
    -n 1 \
    --run_mode "$RUNMODE" \
    --strandedness "$STRAND" \
    --DP "$MINCOV" \
    --AC "$MINALT" \
    --AB "$MINEDIT" \
    --min-base_quality "$QUAL" \
    --min-map_quality "$QUAL" 2>/dev/null || { : > "_snvs_${CHR}.txt"; : > "_bl_${CHR}.txt"; exit 0; }

if [ -s "${PREF}.vcf" ]; then
    : > "_bl_${CHR}.txt"
    scallele_io.py vcf "${PREF}.vcf" "$MINCOV" "_bl_${CHR}.txt" > "_snvs_${CHR}.txt" 2>/dev/null \
        || { : > "_snvs_${CHR}.txt"; : > "_bl_${CHR}.txt"; }
else
    : > "_snvs_${CHR}.txt"; : > "_bl_${CHR}.txt"
fi
rm -f "${PREF}".*
