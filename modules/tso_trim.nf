/*
 * Template-switch oligo (TSO) trimming for paired-end full-length protocols (Smart-seq2), before alignment.
 * A fragment can start in the TSO on either mate, and STAR --clipAdapterType CellRanger4 clips mate 1 only, so the
 * TSO is removed from the 5' end of both mates as a non-anchored 5' adapter (partial TSOs at the read start included).
 * Without it, the TSO's final GGG aligns onto genomic G-A-G and reads as A-to-G.
 */

process CUTADAPT_TSO {
    tag { "${meta.lane_id ?: meta.id}" }
    label 'process_medium'
    conda params.cutadapt_conda

    input:
    tuple val(meta), path(annot), path(r1), path(r2)

    output:
    tuple val(meta), path(annot), path("${meta.lane_id ?: meta.id}_tso_R1.fastq.gz"), path("${meta.lane_id ?: meta.id}_tso_R2.fastq.gz"), emit: reads
    path "${meta.lane_id ?: meta.id}_cutadapt_tso.json", emit: report

    script:
    def prefix = meta.lane_id ?: meta.id
    def tso = params.tso_seq
    """
    cutadapt \\
        -g "${tso};min_overlap=${params.tso_min_overlap}" \\
        -G "${tso};min_overlap=${params.tso_min_overlap}" \\
        --minimum-length ${params.tso_min_length} --pair-filter any \\
        --cores ${task.cpus} \\
        --json ${prefix}_cutadapt_tso.json \\
        -o ${prefix}_tso_R1.fastq.gz -p ${prefix}_tso_R2.fastq.gz \\
        ${r1} ${r2} > ${prefix}_cutadapt_tso.log
    """
}
