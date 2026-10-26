/*
 * GATK preprocessing — SplitNCigarReads and BQSR
 */

process GATK_SPLITNCIGAR {
    tag { "${meta.id}" }
    label 'process_high'
    // Memory escalation: nextflow.config (-Xmx tracks task.memory).
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries    3
    publishDir { "${params.outdir}/${meta.id}/02_gatk_preprocess" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path main_chroms

    output:
    tuple val(meta), path("${meta.id}_split.bam"), path("${meta.id}_split.bam.bai"), emit: bam

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Check if read groups already exist
    HAS_RG=\$(samtools view -H ${bam} | grep -c '^@RG' || true)
    if [ "\$HAS_RG" -eq 0 ]; then
        echo "No read groups found, adding default RG"
        gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}G" AddOrReplaceReadGroups \\
            -I ${bam} \\
            -O _rg.bam \\
            -RGID ${meta.id} \\
            -RGLB lib1 \\
            -RGPL illumina \\
            -RGPU unit1 \\
            -RGSM ${meta.id}
        samtools index _rg.bam
    else
        echo "Read groups present (\$HAS_RG RG lines), keeping existing"
        ln -s ${bam} _rg.bam
        ln -s ${bai} _rg.bam.bai
    fi

    # Ensure FASTA dict and fai exist
    if [ ! -f ${fasta}.fai ]; then samtools faidx ${fasta}; fi
    DICT=\$(echo ${fasta} | sed 's/.fa\$/.dict/' | sed 's/.fasta\$/.dict/')
    if [ ! -f \$DICT ]; then gatk CreateSequenceDictionary -R ${fasta}; fi

    # Restrict traversal to the requested main chromosomes. Externally-aligned
    # BAMs (e.g. Cell Ranger) can carry scaffolds (KI270728.1) with names/order
    # that differ from the analysis-set FASTA dict; SplitNCigarReads iterates
    # reads and aborts on the first read whose contig is absent from the dict.
    LFLAG=""
    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        grep -v '^\$' ${main_chroms} > _intervals.list
        LFLAG="-L _intervals.list"
    fi

    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}G" SplitNCigarReads \\
        -R ${fasta} \\
        -I _rg.bam \\
        \$LFLAG \\
        -O ${meta.id}_split.bam

    samtools index ${meta.id}_split.bam
    rm -f _rg.bam _rg.bam.bai _intervals.list
    """
}


process GATK_BQSR {
    tag { "${meta.id}" }
    label 'process_high'
    // Memory escalation: nextflow.config.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries    3
    publishDir { "${params.outdir}/${meta.id}/02_gatk_preprocess" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path known_sites
    path known_sites_idx
    path main_chroms

    output:
    tuple val(meta), path("${meta.id}_recal.bam"), path("${meta.id}_recal.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_recal_data.table"), emit: table

    when:
    known_sites.name != 'NO_FILE'

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Ensure FASTA dict and fai exist
    if [ ! -f ${fasta}.fai ]; then samtools faidx ${fasta}; fi
    DICT=\$(echo ${fasta} | sed 's/.fa\$/.dict/' | sed 's/.fasta\$/.dict/')
    if [ ! -f \$DICT ]; then gatk CreateSequenceDictionary -R ${fasta}; fi

    # GATK BaseRecalibrator needs a Feature index (.tbi/.idx) beside the known-
    # sites VCF. Nextflow stages the VCF as a symlink, so the index must be
    # staged too (known_sites_idx). Self-heal if it is missing or unstaged.
    if [ ! -f ${known_sites}.tbi ] && [ ! -f ${known_sites}.idx ]; then
        echo "known-sites index not staged; building with IndexFeatureFile" >&2
        gatk IndexFeatureFile -I ${known_sites}
    fi

    # Restrict to main chromosomes (same contig-mismatch guard as SplitNCigar).
    LFLAG=""
    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        grep -v '^\$' ${main_chroms} > _intervals.list
        LFLAG="-L _intervals.list"
    fi

    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}G" BaseRecalibrator \\
        -R ${fasta} \\
        -I ${bam} \\
        \$LFLAG \\
        --known-sites ${known_sites} \\
        -O ${meta.id}_recal_data.table

    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}G" ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        \$LFLAG \\
        --bqsr-recal-file ${meta.id}_recal_data.table \\
        -O ${meta.id}_recal.bam

    samtools index ${meta.id}_recal.bam
    rm -f _intervals.list
    """
}
