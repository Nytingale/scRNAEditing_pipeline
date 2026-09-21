#!/bin/bash
# One-chromosome GATK HaplotypeCaller worker for GATK_HAPLOTYPECALLER.
#
# GATK4 removed HaplotypeCaller's -nct/-nt multithreading, so a whole-genome call
# is a single serial traversal: on the 84 GB Replogle NTC pseudobulk it reached
# only ~31% of the main chromosomes in 22 h (~72 h projected) and would have blown
# the queue wall limit. The module now runs one worker per main chromosome in
# parallel (xargs -P) over the SAME single-read-group BAM and concatenates the
# per-chrom VCFs. This is exact, not an approximation: HaplotypeCaller's local
# reassembly windows never span chromosomes, so a per-contig scatter yields the
# same call set as one traversal.
#
# Args: CHR RGBAM FASTA XMX_GB FILTERS_FILE
# Writes _raw_<CHR>.vcf.gz (+ .tbi) in the current working directory.
set -uo pipefail

CHR="$1"; RGBAM="$2"; FASTA="$3"; XMX_GB="$4"; FILTERS_FILE="$5"
OUT="_raw_${CHR}.vcf.gz"

# The module writes the caller's filter block (which carries backslash-newline
# continuations) to a file; split it back into argv tokens here.
mapfile -t FILTERS < <(tr -s '[:space:]\\' '\n' < "$FILTERS_FILE" | grep -v '^$')

gatk --java-options "-Xmx${XMX_GB}G" HaplotypeCaller \
    -R "$FASTA" \
    -I "$RGBAM" \
    -O "$OUT" \
    -A StrandBiasBySample \
    -L "$CHR" \
    "${FILTERS[@]}"
RC=$?

# Fail loudly: a silently dropped chromosome would quietly shrink the discovered
# site set and corrupt the caller benchmark.
if [ "$RC" -ne 0 ] || [ ! -s "$OUT" ]; then
    echo "ERROR: HaplotypeCaller failed for ${CHR} (exit ${RC})" >&2
    rm -f "$OUT" "${OUT}.tbi"
    exit 1
fi
