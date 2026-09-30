#!/bin/bash
# One-chromosome scAllele worker for SCALLELE (chr-scatter, matches the
# bcftools/gatk pattern). Runs scAllele single-threaded (-n 1) on one contig so
# the module can scatter contigs across the allocated cpus with xargs -P instead
# of the old serial per-chromosome loop. Writes per-chr _snvs_<CHR>.txt (all
# SNVs, PASS + LowQual, with FILTER, varType tag and read-cluster strand; see
# scallele_vcf.py) and keeps the raw per-chr VCF as _vcf_<CHR>.vcf.gz and scAllele's
# read-cluster table as _rc_<CHR>.tsv.gz (index -> strand, span; gives the strand of
# every SNV), so parser changes never need a scAllele rerun. Per-chr output files
# avoid shared-file races.
#
# Args: CHR BAM FASTA MINCOV MINALT MINEDIT QUAL STRAND RUNMODE
# SCALLELE_STRICT=1 (set by SCALLELE_CHROM, one contig per Nextflow task): a
# scAllele or VCF-parse failure exits non-zero so the task fails and retries
# with more memory, instead of silently writing an empty contig.
set -uo pipefail

CHR="$1"; BAM="$2"; FASTA="$3"; MINCOV="$4"; MINALT="$5"; MINEDIT="$6"
QUAL="$7"; STRAND="$8"; RUNMODE="$9"
PREF="_scallele_${CHR}"
STRICT="${SCALLELE_STRICT:-0}"
rm -f "${PREF}".*

empty_out() { : > "_snvs_${CHR}.txt"; : | gzip -c > "_vcf_${CHR}.vcf.gz"; : | gzip -c > "_rc_${CHR}.tsv.gz"; }

fail_or_empty() {
    if [ "$STRICT" = "1" ]; then echo "ERROR: scAllele failed on ${CHR}: $1" >&2; exit 1; fi
    empty_out; exit 0
}

run_scallele() {
    scAllele -b "$BAM" -g "$FASTA" -o "$PREF" -c "$CHR" -n 1 \
        --run_mode "$RUNMODE" --strandedness "$STRAND" --DP "$MINCOV" --AC "$MINALT" \
        --AB "$MINEDIT" --min-base_quality "$QUAL" --min-map_quality "$QUAL"
}

if [ "$STRICT" = "1" ]; then
    run_scallele || fail_or_empty "scAllele exit $?"
else
    run_scallele 2>/dev/null || fail_or_empty "scAllele exit $?"
fi

if [ -s "${PREF}.vcf" ]; then
    RC="${PREF}.read_cluster_info.tab"
    [ -s "$RC" ] || fail_or_empty "no read-cluster table ${RC}"
    scallele_vcf.py "${PREF}.vcf" "$MINCOV" "$RC" > "_snvs_${CHR}.txt" \
        || fail_or_empty "VCF parse exit $?"
    gzip -c "${PREF}.vcf" > "_vcf_${CHR}.vcf.gz"
    gzip -c "$RC" > "_rc_${CHR}.tsv.gz"
else
    empty_out
fi
rm -f "${PREF}".*
exit 0
