/*
 * Strandedness inference using RSeQC
 */

process INFER_STRANDEDNESS {
    tag { "${meta.id}" }
    label 'process_low'
    publishDir { "${params.outdir}/${meta.id}/00_qc" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path refseq_bed

    output:
    tuple val(meta), env('STRANDEDNESS'), emit: strandedness
    tuple val(meta), path("${meta.id}_strandedness.txt"), emit: report

    script:
    """
    #!/bin/bash
    set -euo pipefail

    infer_experiment.py -i ${bam} -r ${refseq_bed} -s 5000000 > _raw.txt 2>/dev/null

    # Parse the two fraction lines. STARsolo 10x/Drop-seq/etc. BAMs are
    # SINGLE-END, so RSeQC prints "++,--" / "+-,-+"; paired-end libraries
    # (e.g. Smart-seq2) print "1++,1--,2+-,2-+" / "1+-,1-+,2++,2--".
    # Match both so each library gets its correct strand.
    FWD=\$(grep -F -e '"1++,1--,2+-,2-+"' -e '"++,--"' _raw.txt | awk '{print \$NF}' | head -1 || true)
    REV=\$(grep -F -e '"1+-,1-+,2++,2--"' -e '"+-,-+"' _raw.txt | awk '{print \$NF}' | head -1 || true)
    [ -z "\$FWD" ] && FWD=0
    [ -z "\$REV" ] && REV=0

    # Determine strandedness
    DIFF=\$(echo "\$FWD - \$REV" | bc)
    ABS_DIFF=\$(echo "\$DIFF" | tr -d '-')

    if (( \$(echo "\$ABS_DIFF < 0.1" | bc -l) )); then
        STRANDEDNESS="unstranded"
        REDITOOLS_S=0
        JACUSA_STRAND="UNSTRANDED"
    elif (( \$(echo "\$DIFF > 0.1" | bc -l) )); then
        STRANDEDNESS="forward"
        REDITOOLS_S=1
        JACUSA_STRAND="FR-SECONDSTRAND"
    else
        STRANDEDNESS="reverse"
        REDITOOLS_S=2
        JACUSA_STRAND="FR-FIRSTSTRAND"
    fi

    echo "=== Strandedness Report: ${meta.id} ===" > ${meta.id}_strandedness.txt
    echo "Forward fraction: \$FWD" >> ${meta.id}_strandedness.txt
    echo "Reverse fraction: \$REV" >> ${meta.id}_strandedness.txt
    echo "Difference: \$DIFF" >> ${meta.id}_strandedness.txt
    echo "Inferred: \$STRANDEDNESS" >> ${meta.id}_strandedness.txt
    echo "REDItools3 -s: \$REDITOOLS_S" >> ${meta.id}_strandedness.txt
    echo "JACUSA2 strand: \$JACUSA_STRAND" >> ${meta.id}_strandedness.txt

    export STRANDEDNESS
    """
}
