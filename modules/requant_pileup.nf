/*
 * Callers compute editing levels differently (QC, BAQ, strand handling), so the
 * published tables get one strand-aware pileup editing level and coverage from the
 * bam_type BAM that produced each call: edit_level_pileup, cov_pileup, edited_reads,
 * ref_reads (original columns kept). Publishes filtered / pre100 to 04_filtered.
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
