/*
 * Collate all downstream filtering results across the benchmarking grid
 * and generate side-by-side comparison plots + summary table.
 */

process COLLATE_RESULTS {
    tag "collate"
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path metrics_tsv

    output:
    path "collated/*", emit: outputs

    script:
    """
    Rscript ${projectDir}/bin/collate_results.R \\
        ${metrics_tsv} \\
        collated
    """
}
