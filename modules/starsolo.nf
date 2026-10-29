/*
 * STARsolo alignment with platform-specific settings (from the benchmarking notebooks).
 */

def getStarsoloParams(platform) {
    def settings = [:]
    settings.soloType = 'CB_UMI_Simple'
    settings.complex  = false
    // 10x 3' TSO; other template-switching platforms use other TSO sequences or paired-end cDNA reads.
    settings.tsoClip  = platform in ['10xv2', '10xv3']

    if (platform == '10xv2') {
        settings.soloCBstart = 1
        settings.soloCBlen   = 16
        settings.soloUMIstart = 17
        settings.soloUMIlen  = 10
        settings.whitelist   = '737K-august-2016.txt'
        settings.hasUMI      = true
    } else if (platform == '10xv3') {
        settings.soloCBstart = 1
        settings.soloCBlen   = 16
        settings.soloUMIstart = 17
        settings.soloUMIlen  = 12
        settings.whitelist   = '3M-february-2018.txt'
        settings.hasUMI      = true
    } else if (platform == 'dropseq') {
        settings.soloCBstart = 1
        settings.soloCBlen   = 12
        settings.soloUMIstart = 13
        settings.soloUMIlen  = 8
        settings.whitelist   = 'None'
        settings.hasUMI      = true
    } else if (platform == 'celseq2') {
        settings.soloCBstart = 7
        settings.soloCBlen   = 6
        settings.soloUMIstart = 1
        settings.soloUMIlen  = 6
        settings.whitelist   = 'None'
        settings.hasUMI      = true
    } else if (platform == 'seqwell') {
        settings.soloCBstart = 1
        settings.soloCBlen   = 12
        settings.soloUMIstart = 13
        settings.soloUMIlen  = 8
        settings.whitelist   = 'None'
        settings.hasUMI      = true
    } else if (platform == 'dronc') {
        settings.soloCBstart = 1
        settings.soloCBlen   = 12
        settings.soloUMIstart = 13
        settings.soloUMIlen  = 8
        settings.whitelist   = 'None'
        settings.hasUMI      = true
    } else if (platform == 'sciseq') {
        settings.soloCBstart = 9
        settings.soloCBlen   = 10
        settings.soloUMIstart = 1
        settings.soloUMIlen  = 8
        settings.whitelist   = 'None'
        settings.hasUMI      = true
    } else if (platform == 'indrop') {
        settings.soloCBstart = 1
        settings.soloCBlen   = 8
        settings.soloUMIstart = 9
        settings.soloUMIlen  = 6
        settings.whitelist   = 'None'
        settings.hasUMI      = true
    } else if (platform == 'indropv2') {
        settings.soloType    = 'CB_UMI_Complex'
        settings.complex     = true
        settings.hasUMI      = true
    } else if (platform == 'smartseq2') {
        settings.soloType    = 'None'
        settings.whitelist   = 'None'
        settings.hasUMI      = false
    } else {
        error "Unknown platform: ${platform}. Supported: 10xv2, 10xv3, dropseq, celseq2, seqwell, indrop, indropv2, sciseq, dronc, smartseq2"
    }

    settings
}


process STAR_GENOMEGENERATE {
    tag "genome"
    label 'process_high'
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    def sjdbOverhang = params.sjdb_overhang ?: 100
    """
    mkdir star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${sjdbOverhang} \\
        --runThreadN ${task.cpus}
    """
}


process STARSOLO_ALIGN {
    tag { "${meta.lane_id ?: meta.id}" }
    label 'process_high'
    publishDir { "${params.outdir}/${meta.id}/01_starsolo" }, mode: 'copy'

    input:
    tuple val(meta), path(celltype_annot), path(r1), path(r2)
    path star_index
    path fasta
    path gtf

    output:
    tuple val(meta), path("${meta.lane_id ?: meta.id}_Aligned.sortedByCoord.out.bam"),
          path("${meta.lane_id ?: meta.id}_Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("${meta.lane_id ?: meta.id}_Solo.out"), emit: solo_out
    path "${meta.lane_id ?: meta.id}_Log.final.out", emit: log

    script:
    def prefix = meta.lane_id ?: meta.id
    def platform = meta.platform
    def settings = getStarsoloParams(platform)

    if (settings.soloType == 'None') {
        // Smart-seq2: regular STAR alignment, no barcode demux
        """
        STAR \\
            --genomeDir ${star_index} \\
            --readFilesIn ${r1} ${r2} \\
            --readFilesCommand zcat \\
            --twopassMode Basic \\
            --outFilterMultimapNmax 1 \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMstrandField intronMotif \\
            --outSAMattributes NH HI nM AS MD \\
            --outFilterType BySJout \\
            --outFileNamePrefix ${prefix}_ \\
            --runThreadN ${task.cpus}

        samtools index ${prefix}_Aligned.sortedByCoord.out.bam
        mkdir -p ${prefix}_Solo.out
        touch ${prefix}_Solo.out/placeholder
        """
    } else if (settings.complex) {
        // InDrop v2: CB_UMI_Complex with 2 gel barcode whitelists
        def default_wl = "${projectDir}/assets/whitelists/inDropCB1.txt ${projectDir}/assets/whitelists/inDropCB2.txt"
        def wl_files = meta.whitelist_files ?: default_wl
        def adapter  = meta.adapter_seq     ?: 'GAGTGATTGCTTGTGACGCCTT'
        def cb_pos   = meta.soloCBposition  ?: '0_0_2_-1 3_1_3_8'
        def umi_pos  = meta.soloUMIposition ?: '3_9_3_14'
        """
        STAR \\
            --genomeDir ${star_index} \\
            --readFilesIn ${r1} ${r2} \\
            --readFilesCommand zcat \\
            --soloType CB_UMI_Complex \\
            --soloCBwhitelist ${wl_files} \\
            --soloCBmatchWLtype 1MM \\
            --soloAdapterSequence ${adapter} \\
            --soloCBposition ${cb_pos} \\
            --soloUMIposition ${umi_pos} \\
            --soloBarcodeReadLength 0 \\
            --soloFeatures GeneFull SJ \\
            --twopassMode Basic \\
            --outFilterMultimapNmax 1 \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN MD \\
            --outFilterType BySJout \\
            --outFileNamePrefix ${prefix}_ \\
            --runThreadN ${task.cpus}

        samtools index ${prefix}_Aligned.sortedByCoord.out.bam
        """
    } else {
        // Standard CB_UMI_Simple platforms
        def whitelist_flag = settings.whitelist != 'None' ?
            "--soloCBwhitelist ${projectDir}/assets/whitelists/${settings.whitelist}" :
            "--soloCBwhitelist None"

        def expected_len = settings.soloCBlen + settings.soloUMIlen
        def clip_flag = (settings.tsoClip && params.tso_clip.toString() != 'false') ?
            "--clipAdapterType CellRanger4" : ""

        """
        #!/bin/bash
        set -euo pipefail

        # Detect actual barcode read length from first record of R1
        if [[ "${r1}" == *.gz ]]; then
            ACTUAL_LEN=\$(nf_head.sh 2 zcat ${r1} | tail -1 | tr -d '\\n' | wc -c)
        else
            ACTUAL_LEN=\$(nf_head.sh 2 cat ${r1} | tail -1 | tr -d '\\n' | wc -c)
        fi

        EXPECTED_LEN=${expected_len}

        if [ "\$ACTUAL_LEN" -ne "\$EXPECTED_LEN" ]; then
            echo "WARNING: [${meta.id}] Platform '${platform}' expects barcode read length of \${EXPECTED_LEN}bp (CB=${settings.soloCBlen} + UMI=${settings.soloUMIlen}) but detected \${ACTUAL_LEN}bp." >&2
            echo "         Using --soloBarcodeReadLength 0 to handle this safely." >&2
        fi

        STAR \\
            --genomeDir ${star_index} \\
            --readFilesIn ${r2} ${r1} \\
            --readFilesCommand zcat \\
            --soloType CB_UMI_Simple \\
            --soloCBstart ${settings.soloCBstart} \\
            --soloCBlen ${settings.soloCBlen} \\
            --soloUMIstart ${settings.soloUMIstart} \\
            --soloUMIlen ${settings.soloUMIlen} \\
            --soloBarcodeReadLength 0 \\
            ${whitelist_flag} \\
            ${clip_flag} \\
            --soloFeatures GeneFull SJ \\
            --twopassMode Basic \\
            --outFilterMultimapNmax 1 \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN MD \\
            --outFilterType BySJout \\
            --outBAMsortingBinsN 200 \\
            --outFileNamePrefix ${prefix}_ \\
            --runThreadN ${task.cpus}

        samtools index ${prefix}_Aligned.sortedByCoord.out.bam
        """
    }
}
