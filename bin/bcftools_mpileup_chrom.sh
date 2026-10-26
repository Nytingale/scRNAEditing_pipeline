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
# Args: CHR BAM FASTA QUAL MINCOV [BAQ]
#   BAQ: nobaq (default) = -B; baq = bcftools default partial BAQ (--bcftools_baq).
# Writes _var_<CHR>.txt and _bl_<CHR>.txt in the current working directory.
set -uo pipefail

CHR="$1"; BAM="$2"; FASTA="$3"; QUAL="$4"; MINCOV="$5"; BAQ="${6:-nobaq}"
PILE="_pile_${CHR}.bcf"
case "$BAQ" in
    nobaq) BAQ_FLAG="-B" ;;
    baq)   BAQ_FLAG="" ;;
    *)     echo "ERROR: BAQ must be baq or nobaq, got '$BAQ'" >&2; exit 1 ;;
esac

# BAMs with more than one @RG SM (Replogle: one per GEM well) must be pooled:
# otherwise mpileup emits one sample per SM and [%DP]/[%AD] below are printed
# per sample concatenated without separator, which downstream_filter.R parses
# as a single garbage DP and the first sample's AD. Single-SM BAMs are unchanged.
N_SM=$(samtools view -H "$BAM" | awk -F'\t' '/^@RG/ {for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) print substr($i, 4)}' | sort -u | wc -l)
IGNORE_RG=""
[ "$N_SM" -gt 1 ] && IGNORE_RG="--ignore-RG"

# Single pileup for this chromosome (-I: skip indels; matches the variant table,
# and the blacklist keeps SNVs only anyway). Written as BCF, read twice below.
bcftools mpileup $IGNORE_RG $BAQ_FLAG \
    -d 100000000 \
    -q "$QUAL" -Q "$QUAL" \
    -r "$CHR" \
    -a FORMAT/AD,FORMAT/DP \
    -I \
    -f "$FASTA" \
    -Ob -o "$PILE" \
    "$BAM" 2> "_pile_${CHR}.err" || { echo "ERROR: mpileup failed for $CHR" >&2; cat "_pile_${CHR}.err" >&2; exit 1; }
rm -f "_pile_${CHR}.err"

# (a) Variant table — identical columns to the original single-pass output.
# Drop non-variant positions (ALT == "<*>" only) at the per-chromosome streaming
# stage — this is the ORIGINAL strategy (see method_scripts/run_mpileup_per_chr.pbs:
# `awk -F'\t' '$5 != "<*>"'`). mpileup emits EVERY covered base with the full
# 16-float I16, so keeping <*>-only rows bloats the whole-arm table to ~9 GB and
# OOMs / breaks fread() in downstream_filter.R. These rows carry no real ALT and
# are dropped anyway at the multi-mismatch step, so this is a lossless reduction
# done here (memory-safe, chr-parallel) instead of loading 9 GB into R.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t%INFO/I16\n' "$PILE" \
    | awk -F'\t' '$4 != "<*>"' > "_var_${CHR}.txt" || { echo "ERROR: bcftools query failed for $CHR" >&2; exit 1; }

# Reads do not map to assembly gaps, so REF=N rows should be ~0. A chromosome that
# is mostly REF=N means mpileup read a broken .fai.
read -r N_ROWS N_REFN < <(awk -F'\t' '{n++} $3 == "N" {m++} END {print n+0, m+0}' "_var_${CHR}.txt")
if [ "$N_ROWS" -ge 100 ] && [ $((N_REFN * 2)) -gt "$N_ROWS" ]; then
    echo "ERROR: $CHR has $N_REFN / $N_ROWS rows with REF=N (broken FASTA index?)" >&2
    exit 1
fi

# (b) Per-sample SNP blacklist: genotype with `call -m`, keep homozygous-ALT
#     (GT="AA") single-base SNVs at DP>=MINCOV. Editing-safe (A-to-I never reaches
#     clean hom-ALT in pseudobulk; het is excluded).
bcftools call -m -Ou "$PILE" | \
    bcftools view -i 'GT="AA"' | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\n' | \
    awk -F'\t' -v mc="$MINCOV" 'length($3)==1 && length($4)==1 && ($5+0)>=mc {print $1"\t"$2}' \
    > "_bl_${CHR}.txt" || { echo "ERROR: SNP-blacklist genotyping failed for $CHR" >&2; exit 1; }

rm -f "$PILE"
