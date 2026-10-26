/*
 * REDIportal overlap: known vs novel sites per (sample, caller, bam_type, strand_mode,
 * snp_sub), all pre100 tables in one R session.
 */

process REDIPORTAL_OVERLAP {
    tag { "rediportal_overlap" }
    label 'process_medium'
    memory { 16.GB * task.attempt }
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path pre100_files, stageAs: 'pre100/*'
    path rediportal

    output:
    path 'rediportal_overlap.tsv', emit: summary

    script:
    """
    Rscript ${projectDir}/bin/rediportal_overlap.R \\
        --rediportal ${rediportal} \\
        --inputs_glob 'pre100/*_pre100.tsv.gz' \\
        --output rediportal_overlap.tsv
    """
}
