/*
 * REQUANT_PILEUP — method-independent re-quantification of editing level and
 * coverage for the published per-site tables.
 *
 * Each caller reports edit_level = cov_alt/(cov_ref+cov_alt) from its own
 * AD/base bookkeeping (different QC/BAQ/strand handling per method), so the
 * values are not strictly comparable and may be trusted verbatim downstream.
 * Here we re-derive one canonical, strand-aware editing level + coverage
 * directly from the SAME bam_type BAM that produced each output, using the
 * project's stranded pileup parser, and ADD columns (originals preserved):
 *     edit_level_pileup, cov_pileup, edited_reads, ref_reads
 *
 * The pileup value is a physical property of the reads, so it is identical
 * across strand_mode variants of a given (sample, bam_type); only which sites
 * PASS differs by strand_mode. It depends on bam_type (dedup/recal change the
 * reads), hence the matching BAM is always used.
 *
 * Publishes the requantified filtered/pre100 under the canonical 04_filtered
 * path (replacing the raw copies, which DOWNSTREAM_R skips when requant is on).
 */

process REQUANT_PILEUP {
    tag { "${meta.id}:${caller}:${bam_label}:${strand_mode}:snp_${snp_sub}" }
    label 'process_low'
    publishDir(
        path: { "${params.outdir}/${meta.id}/04_filtered/${caller}/${bam_label}/${strand_mode}/snp_${snp_sub}" },
        mode: 'copy',
        saveAs: { fn -> fn.startsWith('rq_') ? fn.replaceFirst('^rq_', '') : null }
    )

    input:
    tuple val(meta), val(caller), val(bam_label), val(strand_mode), val(snp_sub),
          path(filtered), path(pre100), path(bam), path(bai), val(libtype)
    path fasta
    path fai

    output:
    tuple val(meta), val(caller), val(bam_label), val(strand_mode), val(snp_sub),
          path("rq_${filtered}"), path("rq_${pre100}"), emit: requant

    script:
    def lt      = (libtype in ['forward', 'reverse', 'unstranded']) ? libtype : 'forward'
    def min_bq  = params.requant_min_bq ?: 30
    def min_mq  = params.requant_min_mq ?: 30
    def maxdep  = params.requant_max_depth ?: 1000000
    """
    set -euo pipefail

    # 1) union of sites across both tables -> 0-based BED
    Rscript ${projectDir}/bin/requant_sites_to_bed.R \\
        --files ${filtered},${pre100} --bed sites.bed

    # 2) single strand-aware mpileup pass at those sites
    if [ -s sites.bed ]; then
        samtools mpileup -B -q ${min_mq} -Q ${min_bq} -d ${maxdep} \\
            -l sites.bed -f ${fasta} ${bam} 2>/dev/null \\
          | perl ${projectDir}/bin/parse_pileup_stranded.pl ${lt} > pileup.tsv
    else
        # no sites: emit just the header so the joiner yields all-NA
        perl ${projectDir}/bin/parse_pileup_stranded.pl ${lt} < /dev/null > pileup.tsv
    fi

    # 3) add canonical columns to each table (originals preserved)
    Rscript ${projectDir}/bin/requant_join.R \\
        --infile ${filtered} --pileup pileup.tsv --outfile rq_${filtered}
    Rscript ${projectDir}/bin/requant_join.R \\
        --infile ${pre100}   --pileup pileup.tsv --outfile rq_${pre100}
    """
}
