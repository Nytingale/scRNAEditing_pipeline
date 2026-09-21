/*
 * BAM validation and repair module
 */

process VALIDATE_BAMS {
    tag { "${meta.lane_id ?: meta.id}" }
    label 'process_low'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.lane_id ?: meta.id}_validation.txt"), emit: report

    script:
    def prefix = meta.lane_id ?: meta.id
    def percell = (meta.platform in ['smartseq2', 'smartseq', 'ss2'])
    """
    #!/bin/bash
    set -euo pipefail

    REPORT="${prefix}_validation.txt"
    PASS=true

    echo "=== BAM Validation Report: ${meta.id} ===" > \$REPORT
    echo "File: ${bam}" >> \$REPORT
    echo "" >> \$REPORT

    # Sample first 10000 reads to check tags
    nf_head.sh 10000 samtools view ${bam} > _sample.sam
    NREADS=\$(wc -l < _sample.sam)
    echo "Sampled \$NREADS reads for tag inspection" >> \$REPORT

    # Check each required tag
    for TAG in CB CR UB MD NH; do
        COUNT=\$(grep -c "\t\${TAG}:" _sample.sam || true)
        PCT=\$(echo "scale=1; \$COUNT * 100 / \$NREADS" | bc 2>/dev/null || echo "0")

        if [ "\$COUNT" -eq 0 ]; then
            STATUS="MISSING"
            if [ "\$TAG" = "MD" ]; then
                echo "  \$TAG: \$STATUS (will be added by samtools calmd)" >> \$REPORT
            elif [ "\$TAG" = "UB" ]; then
                echo "  \$TAG: \$STATUS (WARNING: needed for future UMI-aware dedup)" >> \$REPORT
            elif [ "\$TAG" = "CB" ] || [ "\$TAG" = "CR" ]; then
                if [ "${percell}" = "true" ]; then
                    echo "  \$TAG: \$STATUS (expected for per-cell platform ${meta.platform})" >> \$REPORT
                else
                    echo "  \$TAG: \$STATUS (CRITICAL: required)" >> \$REPORT
                    PASS=false
                fi
            else
                echo "  \$TAG: \$STATUS (CRITICAL: required)" >> \$REPORT
                PASS=false
            fi
        else
            echo "  \$TAG: PRESENT (\${PCT}% of sampled reads)" >> \$REPORT
        fi
    done

    echo "" >> \$REPORT

    # Check if reads are filtered for unique mapping
    MULTI=\$(grep -c "NH:i:[2-9]\\|NH:i:[1-9][0-9]" _sample.sam || true)
    if [ "\$MULTI" -gt 0 ]; then
        PCT_MULTI=\$(echo "scale=1; \$MULTI * 100 / \$NREADS" | bc 2>/dev/null || echo "0")
        echo "Multi-mapped reads: \$MULTI (\${PCT_MULTI}%) — will be filtered" >> \$REPORT
    else
        echo "Multi-mapped reads: none detected (good)" >> \$REPORT
    fi

    echo "" >> \$REPORT

    if [ "\$PASS" = "false" ]; then
        echo "RESULT: FAIL — missing critical tags. Remap with STARsolo (--remap true)." >> \$REPORT
        cat \$REPORT >&2
        exit 1
    else
        echo "RESULT: PASS" >> \$REPORT
    fi

    rm -f _sample.sam
    """
}


process ADD_MD_TAG {
    tag { "${meta.lane_id ?: meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/00_validated" }, mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("${meta.lane_id ?: meta.id}_md.bam"), path("${meta.lane_id ?: meta.id}_md.bam.bai"), emit: bam

    script:
    def prefix = meta.lane_id ?: meta.id
    """
    #!/bin/bash
    set -euo pipefail

    # Check if MD tag already present
    HAS_MD=\$(nf_head.sh 1000 samtools view ${bam} | grep -c "MD:" || true)

    if [ "\$HAS_MD" -gt 0 ]; then
        echo "MD tag already present, creating symlink"
        ln -s ${bam} ${prefix}_md.bam
    else
        echo "MD tag missing, adding with samtools calmd"
        samtools calmd -b ${bam} ${fasta} > ${prefix}_md.bam 2>/dev/null
    fi

    samtools index ${prefix}_md.bam
    """
}


process FILTER_UNIQUELY_MAPPED {
    tag { "${meta.lane_id ?: meta.id}" }
    label 'process_high'
    // Pre-existing BAMs (remap=false) can carry multimappers (e.g. Cell Ranger
    // NH:i:2..N); NOT reached when remap=true (STARsolo already emits unique-only
    // via --outFilterMultimapNmax 1). Large pooled BAMs (~90 GB) previously blew
    // past the 8 h wall as a single-thread awk-over-SAM pass; now a threaded
    // htslib expression filter. Retry on transient SGE kills, escalating time.
    // The time escalation is set in nextflow.config (withName); a `time`
    // directive here is silently overridden by 'process_high'.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries    2
    publishDir { "${params.outdir}/${meta.id}/00_validated" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.lane_id ?: meta.id}_uniq.bam"), path("${meta.lane_id ?: meta.id}_uniq.bam.bai"), emit: bam

    script:
    def prefix = meta.lane_id ?: meta.id
    """
    #!/bin/bash
    set -euo pipefail

    # Keep only uniquely-mapped reads (NH:i:1) via a threaded htslib expression
    # filter — orders of magnitude faster than parsing SAM text with awk. Input
    # is coordinate-sorted (from ADD_MD_TAG) and the filter preserves order, so
    # no re-sort is needed; just index the output.
    samtools view -b -e '[NH]==1' -@ ${task.cpus} -o ${prefix}_uniq.bam ${bam}
    samtools index -@ ${task.cpus} ${prefix}_uniq.bam

    # Read counts come straight from the .bai (idxstats), avoiding full extra
    # passes over the (large) BAMs.
    BEFORE=\$(samtools idxstats ${bam} | awk '{s+=\$3+\$4} END{print s}')
    AFTER=\$(samtools idxstats ${prefix}_uniq.bam | awk '{s+=\$3+\$4} END{print s}')
    PCT=\$(awk -v a=\$AFTER -v b=\$BEFORE 'BEGIN{ if(b>0) printf "%.1f", a*100/b; else print "NA" }')
    echo "Unique mapping filter: \$BEFORE -> \$AFTER reads (\${PCT}%)"
    """
}
