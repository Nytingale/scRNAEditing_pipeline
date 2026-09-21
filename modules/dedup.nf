/*
 * Deduplication module — three strategies:
 *   'none'       — skip dedup entirely
 *   'coordinate' — cell-aware coordinate-based (GATK MarkDuplicates)
 *   'umi'        — UMI-aware (umi_tools dedup)
 *
 * Cell barcode tag is configurable via params.dedup_cell_tag (default: CR).
 * CR (raw barcode) is recommended because CB can be "-" for uncorrectable
 * barcodes, which would lump unrelated reads into one group.
 * Users who prefer CB can set --dedup_cell_tag CB; a guard step automatically
 * removes reads with CB:Z:- (only when using CB, to avoid unnecessary I/O).
 */

process GATK_MARKDUPLICATES {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/02_dedup" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_rmDup.bam"), path("${meta.id}_rmDup.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_dup_metrics.txt"), emit: metrics

    script:
    def cell_tag = params.dedup_cell_tag ?: 'CR'
    """
    #!/bin/bash
    set -euo pipefail

    INPUT_BAM="${bam}"

    # When using CB, guard against CB:Z:- (uncorrectable barcodes).
    # Skipped for CR since raw barcodes are always real sequences.
    if [ "${cell_tag}" = "CB" ]; then
        N_TOTAL=\$(samtools view -c ${bam})
        samtools view -h ${bam} | awk '/^@/ || \$0 !~ /CB:Z:-/' | \\
            samtools view -bS -o _clean.bam -
        samtools index _clean.bam
        N_CLEAN=\$(samtools view -c _clean.bam)
        N_DROPPED=\$((N_TOTAL - N_CLEAN))
        if [ "\$N_DROPPED" -gt 0 ]; then
            echo "Removed \${N_DROPPED} reads with CB:Z:- before dedup (\${N_CLEAN} remain)"
        fi
        INPUT_BAM="_clean.bam"
    fi

    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}G" MarkDuplicates \\
        -I \${INPUT_BAM} \\
        -O ${meta.id}_rmDup.bam \\
        -M ${meta.id}_dup_metrics.txt \\
        --BARCODE_TAG ${cell_tag} \\
        --REMOVE_DUPLICATES true

    samtools index ${meta.id}_rmDup.bam
    rm -f _clean.bam _clean.bam.bai
    """
}


process UMI_DEDUP {
    tag { "${meta.id}" }
    label 'process_high'
    publishDir { "${params.outdir}/${meta.id}/02_dedup" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_umiDup.bam"), path("${meta.id}_umiDup.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_umi_dedup_stats*"), emit: stats

    script:
    def cell_tag = params.dedup_cell_tag ?: 'CR'
    """
    #!/bin/bash
    set -euo pipefail

    # Validate UB tag presence
    HAS_UB=\$(nf_head.sh 1000 samtools view ${bam} | grep -c "UB:Z:" || true)
    if [ "\$HAS_UB" -eq 0 ]; then
        echo "ERROR: UMI-aware dedup requires UB (corrected UMI) tag but none found." >&2
        echo "This BAM may be from a Smart-seq2 or plate-based protocol without UMIs." >&2
        echo "Use --dedup coordinate or --dedup none instead." >&2
        exit 1
    fi

    INPUT_BAM="${bam}"

    # When using CB, guard against CB:Z:- (uncorrectable barcodes).
    # Skipped for CR since raw barcodes are always real sequences.
    if [ "${cell_tag}" = "CB" ]; then
        N_TOTAL=\$(samtools view -c ${bam})
        samtools view -h ${bam} | awk '/^@/ || \$0 !~ /CB:Z:-/' | \\
            samtools view -bS -o _clean.bam -
        samtools index _clean.bam
        N_CLEAN=\$(samtools view -c _clean.bam)
        N_DROPPED=\$((N_TOTAL - N_CLEAN))
        if [ "\$N_DROPPED" -gt 0 ]; then
            echo "Removed \${N_DROPPED} reads with CB:Z:- before UMI dedup (\${N_CLEAN} remain)"
        fi
        INPUT_BAM="_clean.bam"
    fi

    umi_tools dedup \\
        --per-cell \\
        --cell-tag=${cell_tag} \\
        --umi-tag=UB \\
        --extract-umi-method=tag \\
        --method=unique \\
        --output-stats=${meta.id}_umi_dedup_stats \\
        -I \${INPUT_BAM} \\
        -S ${meta.id}_umiDup.bam \\
        --log=${meta.id}_umi_dedup.log

    samtools index ${meta.id}_umiDup.bam
    rm -f _clean.bam _clean.bam.bai

    # Report dedup stats
    BEFORE=\$(samtools view -c ${bam})
    AFTER=\$(samtools view -c ${meta.id}_umiDup.bam)
    echo "UMI dedup: \${BEFORE} -> \${AFTER} reads (\$(echo "scale=1; \$AFTER * 100 / \$BEFORE" | bc)% retained)"
    """
}


/*
 * Passthrough process for --dedup none
 * Simply symlinks the input BAM with a consistent naming scheme
 */
process SKIP_DEDUP {
    tag { "${meta.id}" }
    label 'process_low'
    publishDir { "${params.outdir}/${meta.id}/02_dedup" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_noDup.bam"), path("${meta.id}_noDup.bam.bai"), emit: bam

    script:
    """
    cp -L ${bam} ${meta.id}_noDup.bam
    cp -L ${bai} ${meta.id}_noDup.bam.bai
    """
}
