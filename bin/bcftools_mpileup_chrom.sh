#!/bin/bash
# One-chromosome bcftools pileup worker for BCFTOOLS_MPILEUP.
#
# A SINGLE mpileup pass per chromosome now feeds BOTH outputs (previously two
# separate full-genome pileups — one for the variant table, one for the SNP
# blacklist). The pileup (the expensive step) is written once to a per-chrom
# BCF, then reused for (a) the raw variant table and (b) hom-ALT SNV genotyping.
# Called in parallel across chromosomes (xargs -P) by the module, so the old
# serial per-chrom loop + single-threaded mpileup no longer serializes the run.
#
# Args: CHR BAM FASTA QUAL MINCOV
# Writes _var_<CHR>.txt and _bl_<CHR>.txt in the current working directory.
set -uo pipefail

CHR="$1"; BAM="$2"; FASTA="$3"; QUAL="$4"; MINCOV="$5"
PILE="_pile_${CHR}.bcf"

# Single pileup for this chromosome (-I: skip indels; matches the variant table,
# and the blacklist keeps SNVs only anyway). Written as BCF, read twice below.
bcftools mpileup \
    -d 100000000 \
    -q "$QUAL" -Q "$QUAL" \
    -r "$CHR" \
    -a FORMAT/AD,FORMAT/DP \
    -I \
    -f "$FASTA" \
    -Ob -o "$PILE" \
    "$BAM" 2>/dev/null || { echo "WARN: mpileup failed for $CHR" >&2; : > "_var_${CHR}.txt"; : > "_bl_${CHR}.txt"; exit 0; }

# (a) Variant table — identical columns to the original single-pass output.
# Drop non-variant positions (ALT == "<*>" only) at the per-chromosome streaming
# stage — this is the ORIGINAL strategy (see method_scripts/run_mpileup_per_chr.pbs:
# `awk -F'\t' '$5 != "<*>"'`). mpileup emits EVERY covered base with the full
# 16-float I16, so keeping <*>-only rows bloats the whole-arm table to ~9 GB and
# OOMs / breaks fread() in downstream_filter.R. These rows carry no real ALT and
# are dropped anyway at the multi-mismatch step, so this is a lossless reduction
# done here (memory-safe, chr-parallel) instead of loading 9 GB into R.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t%INFO/I16\n' "$PILE" 2>/dev/null \
    | awk -F'\t' '$4 != "<*>"' > "_var_${CHR}.txt" || : > "_var_${CHR}.txt"

# (b) Per-sample SNP blacklist: genotype with `call -m`, keep homozygous-ALT
#     (GT="AA") single-base SNVs at DP>=MINCOV. Editing-safe (A-to-I never reaches
#     clean hom-ALT in pseudobulk; het is excluded).
bcftools call -m -Ou "$PILE" 2>/dev/null | \
    bcftools view -i 'GT="AA"' 2>/dev/null | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\n' 2>/dev/null | \
    awk -F'\t' -v mc="$MINCOV" 'length($3)==1 && length($4)==1 && ($5+0)>=mc {print $1"\t"$2}' \
    > "_bl_${CHR}.txt" 2>/dev/null || : > "_bl_${CHR}.txt"

rm -f "$PILE"
