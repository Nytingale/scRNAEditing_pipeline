/*
 * Downstream filtering: gene annotation, SNP removal, strand assignment and removal
 * of sites near other mismatch types.
 */

// `--flag false` arrives as the string "false", which is truthy in Groovy.
def toBool(v) {
    if (v instanceof Boolean) return v
    return (v?.toString()?.trim()?.toLowerCase() in ['true', '1', 'yes', 'y', 't', 'on'])
}

process DOWNSTREAM_R {
    tag { "${meta.id}:${caller}:${bam_label}:${strand_mode}:snp_${snp_sub}" }
    label 'process_medium'
    // Memory ladder: nextflow.config. Exit 3 = scAllele 'caller' strand mode
    // unavailable (unstranded library): the configuration does not exist, so no retry.
    errorStrategy { task.exitStatus == 3 ? 'ignore' : (task.attempt <= 3 ? 'retry' : 'finish') }
    maxRetries    3
    publishDir(
        path: { "${params.outdir}/${meta.id}/04_filtered/${caller}/${bam_label}/${strand_mode}/snp_${snp_sub}" },
        mode: 'copy',
        // With requant on, REQUANT_PILEUP publishes filtered / pre100.
        saveAs: { fn -> (toBool(params.requant_editing) && (fn.endsWith('_filtered.txt.gz') || fn.endsWith('_pre100.tsv.gz'))) ? null : fn }
    )

    input:
    tuple val(meta), val(caller), val(bam_label), path(variants), val(strand_mode), path(snp_blacklist), val(snp_sub)
    path fasta
    path gtf
    path snp_file
    path strand_annots, stageAs: 'strand_annots/*'
    path alu_annots,    stageAs: 'alu_annots/*'
    path main_chroms

    output:
    tuple val(meta), val(caller), val(bam_label), val(strand_mode), val(snp_sub),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_filtered.txt.gz"), emit: filtered
    tuple val(meta), val(caller), val(bam_label), val(strand_mode), val(snp_sub),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_pre100.tsv.gz"), emit: pre100
    tuple val(meta), val(caller), val(bam_label), val(strand_mode),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_summary.txt"), emit: summary
    tuple val(meta), val(caller), val(bam_label), val(strand_mode),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_stepwise.tsv"), emit: stepwise
    tuple val(meta), val(caller), val(bam_label), val(strand_mode),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_metrics.tsv"), emit: metrics

    script:
    def snp_flag = ""
    if (snp_file.name != 'NO_FILE') {
        def ext = snp_file.name.toLowerCase()
        if (ext.endsWith('.rdata') || ext.endsWith('.rds')) {
            snp_flag = "--snp_rdata ${snp_file}"
        } else {
            snp_flag = "--snp_bed ${snp_file}"
        }
    }
    def blacklist_flag = snp_blacklist.name != 'NO_FILE' ? "--snp_blacklist ${snp_blacklist}" : ""
    def strand_flag  = strand_annots.name != 'NO_FILE' ? "--strand_annots ${strand_annots}" : ""
    def alu_flag     = alu_annots.name   != 'NO_FILE' ? "--alu_annots ${alu_annots}" : ""
    def chroms_flag  = main_chroms.name  != 'NO_FILE' ? "--main_chroms ${main_chroms}" : ""
    def windows      = params.adaptive_windows
    def mm_window    = toBool(params.mismatch_window_filter) ? 1 : 0
    def clustered    = toBool(params.clustered_filter) ? 1 : 0
    def save_int_n   = (params.save_intermediates ?: 0).toString().isInteger() ? (params.save_intermediates).toString().toInteger() : 0
    def intermed     = save_int_n > 0 ? "--save_intermediates ${save_int_n}" : ""
    def ambig_thresh = params.strand_ambiguity_threshold
    def min_alt      = params.min_alt
    def min_edit     = params.min_edit_level
    // 'count': cov_alt >= min_alt; 'binom': alt count above the sequencing error rate.
    def alt_rule     = params.alt_rule
    def seq_err      = params.seq_error_rate
    def alt_alpha    = params.alt_binom_alpha
    def alt_floor    = params.min_alt_floor
    def j2_min_score = params.jacusa2_min_score
    def rm_100pct    = toBool(params.remove_100pct_editing) ? 'TRUE' : 'FALSE'
    """
    # downstream_filter v2 — requires R.utils for gzipped I/O; requant pipeline
    Rscript ${projectDir}/bin/downstream_filter_v8.R \\
        --input ${variants} \\
        --caller ${caller} \\
        --sample ${meta.id} \\
        --bam_label ${bam_label} \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf} \\
        ${snp_flag} \\
        ${blacklist_flag} \\
        ${strand_flag} \\
        ${alu_flag} \\
        ${chroms_flag} \\
        --windows "${windows}" \\
        --mismatch_window_filter ${mm_window} \\
        --clustered_filter ${clustered} \\
        --min_cov ${params.min_cov} \\
        --min_alt ${min_alt} \\
        --min_edit_level ${min_edit} \\
        --alt_rule ${alt_rule} \\
        --seq_error_rate ${seq_err} \\
        --alt_binom_alpha ${alt_alpha} \\
        --min_alt_floor ${alt_floor} \\
        --jacusa2_min_score ${j2_min_score} \\
        --remove_100pct_editing ${rm_100pct} \\
        --strand_mode ${strand_mode} \\
        --strand_ambiguity_threshold ${ambig_thresh} \\
        --snp_subtract_label ${snp_sub} \\
        --output_filtered ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_filtered.txt.gz \\
        --output_pre100 ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_pre100.tsv.gz \\
        --output_summary ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_summary.txt \\
        --output_stepwise ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_stepwise.tsv \\
        --output_metrics ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_metrics.tsv \\
        ${intermed}
    """
}

/*
 * Same as DOWNSTREAM_R but runs bin/downstream_filter_scallele.R. A separate process,
 * so routing scAllele does not change the task hash of the other callers.
 */
process DOWNSTREAM_R_SCALLELE {
    tag { "${meta.id}:${caller}:${bam_label}:${strand_mode}:snp_${snp_sub}" }
    label 'process_medium'
    // Memory ladder: nextflow.config. Exit 3 = scAllele 'caller' strand mode
    // unavailable (unstranded library): the configuration does not exist, so no retry.
    errorStrategy { task.exitStatus == 3 ? 'ignore' : (task.attempt <= 3 ? 'retry' : 'finish') }
    maxRetries    3
    publishDir(
        path: { "${params.outdir}/${meta.id}/04_filtered/${caller}/${bam_label}/${strand_mode}/snp_${snp_sub}" },
        mode: 'copy',
        // With requant on, REQUANT_PILEUP publishes filtered / pre100.
        saveAs: { fn -> (toBool(params.requant_editing) && (fn.endsWith('_filtered.txt.gz') || fn.endsWith('_pre100.tsv.gz'))) ? null : fn }
    )

    input:
    tuple val(meta), val(caller), val(bam_label), path(variants), val(strand_mode), path(snp_blacklist), val(snp_sub)
    path fasta
    path gtf
    path snp_file
    path strand_annots, stageAs: 'strand_annots/*'
    path alu_annots,    stageAs: 'alu_annots/*'
    path main_chroms

    output:
    tuple val(meta), val(caller), val(bam_label), val(strand_mode), val(snp_sub),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_filtered.txt.gz"), emit: filtered
    tuple val(meta), val(caller), val(bam_label), val(strand_mode), val(snp_sub),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_pre100.tsv.gz"), emit: pre100
    tuple val(meta), val(caller), val(bam_label), val(strand_mode),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_summary.txt"), emit: summary
    tuple val(meta), val(caller), val(bam_label), val(strand_mode),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_stepwise.tsv"), emit: stepwise
    tuple val(meta), val(caller), val(bam_label), val(strand_mode),
          path("${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_metrics.tsv"), emit: metrics

    script:
    def snp_flag = ""
    if (snp_file.name != 'NO_FILE') {
        def ext = snp_file.name.toLowerCase()
        if (ext.endsWith('.rdata') || ext.endsWith('.rds')) {
            snp_flag = "--snp_rdata ${snp_file}"
        } else {
            snp_flag = "--snp_bed ${snp_file}"
        }
    }
    def blacklist_flag = snp_blacklist.name != 'NO_FILE' ? "--snp_blacklist ${snp_blacklist}" : ""
    def strand_flag  = strand_annots.name != 'NO_FILE' ? "--strand_annots ${strand_annots}" : ""
    def alu_flag     = alu_annots.name   != 'NO_FILE' ? "--alu_annots ${alu_annots}" : ""
    def chroms_flag  = main_chroms.name  != 'NO_FILE' ? "--main_chroms ${main_chroms}" : ""
    def windows      = params.adaptive_windows
    def mm_window    = toBool(params.mismatch_window_filter) ? 1 : 0
    def clustered    = toBool(params.clustered_filter) ? 1 : 0
    def save_int_n   = (params.save_intermediates ?: 0).toString().isInteger() ? (params.save_intermediates).toString().toInteger() : 0
    def intermed     = save_int_n > 0 ? "--save_intermediates ${save_int_n}" : ""
    def ambig_thresh = params.strand_ambiguity_threshold
    def min_alt      = params.min_alt
    def min_edit     = params.min_edit_level
    // 'count': cov_alt >= min_alt; 'binom': alt count above the sequencing error rate.
    def alt_rule     = params.alt_rule
    def seq_err      = params.seq_error_rate
    def alt_alpha    = params.alt_binom_alpha
    def alt_floor    = params.min_alt_floor
    def j2_min_score = params.jacusa2_min_score
    def rm_100pct    = toBool(params.remove_100pct_editing) ? 'TRUE' : 'FALSE'
    """
    # downstream_filter v2 — requires R.utils for gzipped I/O; requant pipeline
    Rscript ${projectDir}/bin/downstream_filter_scallele.R \\
        --input ${variants} \\
        --caller ${caller} \\
        --sample ${meta.id} \\
        --bam_label ${bam_label} \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf} \\
        ${snp_flag} \\
        ${blacklist_flag} \\
        ${strand_flag} \\
        ${alu_flag} \\
        ${chroms_flag} \\
        --windows "${windows}" \\
        --mismatch_window_filter ${mm_window} \\
        --clustered_filter ${clustered} \\
        --min_cov ${params.min_cov} \\
        --min_alt ${min_alt} \\
        --min_edit_level ${min_edit} \\
        --alt_rule ${alt_rule} \\
        --seq_error_rate ${seq_err} \\
        --alt_binom_alpha ${alt_alpha} \\
        --min_alt_floor ${alt_floor} \\
        --jacusa2_min_score ${j2_min_score} \\
        --remove_100pct_editing ${rm_100pct} \\
        --strand_mode ${strand_mode} \\
        --strand_ambiguity_threshold ${ambig_thresh} \\
        --snp_subtract_label ${snp_sub} \\
        --output_filtered ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_filtered.txt.gz \\
        --output_pre100 ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_pre100.tsv.gz \\
        --output_summary ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_summary.txt \\
        --output_stepwise ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_stepwise.tsv \\
        --output_metrics ${meta.id}_${caller}_${bam_label}_${strand_mode}_snp${snp_sub}_metrics.tsv \\
        ${intermed}
    """
}
