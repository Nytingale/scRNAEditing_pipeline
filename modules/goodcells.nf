/*
 * goodCells construction
 *
 * goodCells = annotated ∪ QC-passing unannotated
 *
 * Candidate pool is the union of STARsolo/CellRanger filtered barcodes and
 * annotated barcodes (prevents hierarchy inversion when the caller is more
 * conservative than the annotation). Annotated cells are assumed QC-passed.
 * Unannotated cells are filtered with uniform nFeature / nCount / pct_mt
 * thresholds. Mitochondrial genes are taken from the GTF by chromosome.
 */

process COMPUTE_GOOD_CELLS {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/02_goodcells" }, mode: 'copy'

    input:
    tuple val(meta), path(annot), val(called_path), val(qc_path), val(matrix_dirs), val(cell_ids)
    path gtf

    output:
    tuple val(meta), path("${meta.id}_goodcells_barcodes.tsv"), emit: barcodes
    tuple val(meta), path("${meta.id}_barcode_qc_summary.tsv"), emit: summary
    tuple val(meta), path("${meta.id}_unannotated_qc_metrics.tsv"), emit: metrics

    script:
    def called_arg  = called_path  ? "--called ${called_path}" : ""
    def qc_arg      = qc_path      ? "--qc_barcodes ${qc_path}" : ""
    def matrix_arg  = matrix_dirs  ? "--matrix_dir ${matrix_dirs}" : ""
    def mito_arg    = params.mito_chr ? "--mito_chr ${params.mito_chr}" : ""
    def cell_list   = cell_ids instanceof List ? cell_ids.findAll { it } : (cell_ids ? [cell_ids] : [])
    def cells_arg   = cell_list ? "--sample_cells _sample_cells.txt" : ""
    def cells_block = cell_list ? "printf '%s\\n' ${cell_list.collect { "'${it}'" }.join(' ')} > _sample_cells.txt" : "true"
    """
    #!/bin/bash
    set -euo pipefail

    ${cells_block}

    compute_goodcells.R \\
        --sample ${meta.id} \\
        --annot ${annot} \\
        --gtf ${gtf} \\
        ${called_arg} \\
        ${qc_arg} \\
        ${matrix_arg} \\
        ${cells_arg} \\
        ${mito_arg} \\
        --qc_min_features ${params.qc_min_features} \\
        --qc_max_features ${params.qc_max_features} \\
        --qc_max_pct_mt ${params.qc_max_pct_mt} \\
        --qc_min_counts ${params.qc_min_counts} \\
        --output_barcodes ${meta.id}_goodcells_barcodes.tsv \\
        --output_summary ${meta.id}_barcode_qc_summary.tsv \\
        --output_metrics ${meta.id}_unannotated_qc_metrics.tsv
    """
}
