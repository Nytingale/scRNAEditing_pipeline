/*
 * Deduplication: 'none', 'coordinate' (GATK MarkDuplicates, cell-aware) or 'umi'
 * (umi_tools). CR is the default cell tag: CB can be "-" for uncorrectable barcodes,
 * which would lump unrelated reads; with CB, reads with CB:Z:- are removed first.
 */

/*
 * Shared driver for both strategies: per region (one main chromosome, `_other` or
 * `_all`) and per GEM well. A cell barcode is unique only within a GEM well, so a BAM
 * pooling wells (>1 @RG SM) is split by well, each well deduplicated on its own, then
 * merged. The caller defines `dedup_unit IN.bam OUT.bam` before this block.
 */
def dedup_per_well_script(bam, chrom, main_chroms, cell_tag, cpus, out) {
    """
    if [ "${chrom}" = "_other" ]; then
        samtools idxstats ${bam} | awk -F'\\t' 'NR==FNR { m[\$1] = 1; next } \$1 != "*" && \$3 > 0 && !(\$1 in m) { print \$1 }' ${main_chroms} - > _regions.txt
    elif [ "${chrom}" = "_all" ]; then
        samtools idxstats ${bam} | awk -F'\\t' '\$1 != "*" && \$3 > 0 { print \$1 }' > _regions.txt
    else
        echo "${chrom}" > _regions.txt
    fi
    N_IN=0; N_WELLS=0
    if [ -s _regions.txt ]; then
        samtools view -b -@ ${cpus} -o _in.bam ${bam} \$(tr '\\n' ' ' < _regions.txt)
        N_IN=\$(samtools view -c _in.bam)
    fi
    if [ "\$N_IN" -eq 0 ]; then
        samtools view -H -b -o ${out} ${bam}
    else
        if [ "${cell_tag}" = "CB" ]; then   # CB:Z:- = uncorrectable barcode, never a cell
            samtools view -h _in.bam | awk '/^@/ || \$0 !~ /CB:Z:-/' | samtools view -b -o _clean.bam -
            mv _clean.bam _in.bam
        fi
        if [ "${out.endsWith('_umiDup.bam')}" = "true" ]; then   # UB:Z:- = uncorrectable UMI (STARsolo): cannot be UMI-deduplicated
            samtools view -h _in.bam | awk '/^@/ || \$0 !~ /\\tUB:Z:-(\\t|\$)/' | samtools view -b -o _clean.bam -
            mv _clean.bam _in.bam
        fi
        rm -rf _units _dd; mkdir -p _units _dd
        # RG ID -> SM (well), from the header
        samtools view -H _in.bam | awk -F'\\t' '/^@RG/ { id = ""; sm = ""; for (i = 2; i <= NF; i++) { if (\$i ~ /^ID:/) id = substr(\$i, 4); if (\$i ~ /^SM:/) sm = substr(\$i, 4) }; print id "\\t" sm }' > _rg2sm.tsv
        N_SM=\$(cut -f2 _rg2sm.tsv | sort -u | wc -l)
        if [ "\$N_SM" -gt 1 ]; then
            samtools split -@ ${cpus} -u _units/_noRG.bam -f '_units/rg_%#.bam' _in.bam
            if [ -s _units/_noRG.bam ] && [ "\$(samtools view -c _units/_noRG.bam)" -gt 0 ]; then
                echo "ERROR: reads without @RG in a multi-well BAM (${chrom}); cannot assign a well" >&2; exit 1
            fi
            rm -f _units/_noRG.bam
            # one file per WELL: join the read groups of each SM
            for f in _units/rg_*.bam; do
                id=\$(samtools view -H "\$f" | awk -F'\\t' '/^@RG/ { for (i = 2; i <= NF; i++) if (\$i ~ /^ID:/) print substr(\$i, 4) }' | head -1)
                n=\$(samtools view -c "\$f")
                [ "\$n" -gt 0 ] || { rm -f "\$f"; continue; }
                sm=\$(awk -F'\\t' -v i="\$id" '\$1 == i { print \$2 }' _rg2sm.tsv | head -1 | tr -c 'A-Za-z0-9._\\n-' '_')
                mkdir -p "_units/well_\$sm"; mv "\$f" "_units/well_\$sm/"
            done
            for d in _units/well_*; do
                w=\$(basename "\$d")
                if [ "\$(ls "\$d"/*.bam | wc -l)" -eq 1 ]; then mv "\$d"/*.bam "_units/\$w.bam"
                else samtools merge -c -p -f -o "_units/\$w.bam" "\$d"/*.bam; fi
                rm -rf "\$d"
            done
        else
            mv _in.bam _units/well_all.bam
        fi
        rm -f _in.bam _in.bam.bai
        for u in _units/*.bam; do samtools index "\$u"; done
        N_WELLS=\$(ls _units/*.bam | wc -l)
        export -f dedup_unit
        ls _units/*.bam | xargs -P ${cpus} -I{} bash -c 'set -euo pipefail; dedup_unit "\$1" "_dd/\$(basename "\$1")"' _ {}
        if [ "\$N_WELLS" -eq 1 ]; then mv _dd/*.bam ${out}
        else samtools merge -@ ${cpus} -c -p -f -o ${out} _dd/*.bam; fi
        rm -rf _units
    fi
    echo "dedup ${chrom}: \${N_IN} reads, \${N_WELLS} well(s) deduplicated separately -> \$(samtools view -c ${out}) reads"
    """
}

/*
 * Coordinate dedup per chromosome and GEM well. 10x reads are single-end, so no
 * duplicate set spans two chromosomes.
 */
process MARKDUP_CHROM {
    tag { "${meta.id}:${chrom}" }
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), val(chrom)
    path main_chroms

    output:
    tuple val(meta), val(chrom), path("${meta.id}__${chrom}_rmDup.bam"), path("${meta.id}__${chrom}_dup_metrics.txt"), emit: bam

    script:
    def cell_tag = params.dedup_cell_tag ?: 'CR'
    def out = "${meta.id}__${chrom}_rmDup.bam"
    def met = "${meta.id}__${chrom}_dup_metrics.txt"
    // wells run task.cpus-way in parallel, one JVM each: split the heap
    def heap_gb = Math.max(1, ((task.memory.toGiga() * 0.8) / task.cpus).intValue())
    """
    #!/bin/bash
    set -euo pipefail
    dedup_unit() {
        gatk --java-options "-Xmx${heap_gb}G" MarkDuplicates \\
            -I "\$1" -O "\$2" -M "\${2%.bam}.metrics.txt" \\
            --BARCODE_TAG ${cell_tag} --REMOVE_DUPLICATES true \\
            --QUIET true --VERBOSITY WARNING
    }
    ${dedup_per_well_script(bam, chrom, main_chroms, cell_tag, task.cpus, out)}
    # one Picard section per well; MARKDUP_GATHER sums them per library
    cat _dd/*.metrics.txt > ${met} 2>/dev/null || : > ${met}
    rm -rf _dd
    """
}

// Picard metrics are summed per LIBRARY across chromosomes; ESTIMATED_LIBRARY_SIZE
// is left blank (not additive).
process MARKDUP_GATHER {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/02_dedup" }, mode: 'copy'

    input:
    tuple val(meta), path(bams), path(metrics)

    output:
    tuple val(meta), path("${meta.id}_rmDup.bam"), path("${meta.id}_rmDup.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_dup_metrics.txt"), emit: metrics

    script:
    def bl = (bams instanceof List ? bams : [bams]).collect { it.name }.sort().join(' ')
    def ml = (metrics instanceof List ? metrics : [metrics]).collect { it.name }.sort().join(' ')
    """
    #!/bin/bash
    set -euo pipefail
    SUM=0
    for b in ${bl}; do SUM=\$((SUM + \$(samtools view -c \$b))); done
    samtools merge -@ ${task.cpus} -c -p -f -o ${meta.id}_rmDup.bam ${bl}
    samtools index -@ ${task.cpus} ${meta.id}_rmDup.bam
    N=\$(samtools view -c ${meta.id}_rmDup.bam)
    [ "\$N" -eq "\$SUM" ] || { echo "ERROR: merged \$N reads != sum of pieces \$SUM" >&2; exit 1; }

    awk -F'\\t' -v OFS='\\t' '
        /^## METRICS CLASS/ { inm = 1; if (!hdr_done) { cls = \$0 }; getline; if (!hdr_done) { cols = \$0; nc = split(\$0, h, "\\t") }; next }
        inm && NF == 0 { inm = 0; hdr_done = 1; next }
        inm && \$1 != "" { lib = \$1; if (!(lib in seen)) { seen[lib] = 1; libs[++nl] = lib }
                          for (i = 2; i <= nc; i++) if (h[i] ~ /(READS|PAIRS|RDS|DUPLICATES)/) s[lib, i] += \$i; next }
        END {
            print "## htsjdk.samtools.metrics.StringHeader"
            print "# Summed over per-chromosome MarkDuplicates runs (MARKDUP_GATHER); ESTIMATED_LIBRARY_SIZE not additive, left blank"
            print ""; print cls; print cols
            for (l = 1; l <= nl; l++) {
                lib = libs[l]; line = lib
                for (i = 2; i <= nc; i++) {
                    v = ""
                    if (h[i] ~ /(READS|PAIRS|RDS|DUPLICATES)/) v = s[lib, i] + 0
                    else if (h[i] == "PERCENT_DUPLICATION") { ex = s[lib, ci("UNPAIRED_READS_EXAMINED")] + 2 * s[lib, ci("READ_PAIRS_EXAMINED")]
                        du = s[lib, ci("UNPAIRED_READ_DUPLICATES")] + 2 * s[lib, ci("READ_PAIR_DUPLICATES")]
                        v = (ex > 0) ? sprintf("%.6f", du / ex) : 0 }
                    line = line OFS v
                }
                print line
            }
            print ""
        }
        function ci(n,   k) { for (k = 1; k <= nc; k++) if (h[k] == n) return k; return 0 }
    ' ${ml} > ${meta.id}_dup_metrics.txt
    echo "MarkDuplicates ${meta.id}: \$N reads from \$(echo ${bl} | wc -w) chromosome pieces"
    """
}


/*
 * Smart-seq2: reads carry no cell tag, so each per-cell BAM is deduplicated on its own
 * (a merged-BAM dedup would collapse identical fragments of different cells).
 * Outputs are named by zero-padded input position so they stay parallel to cell_ids.
 */
process MARKDUP_PERCELL_SS2 {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/02_dedup" }, mode: 'copy', pattern: '*_{rmDup.bam,rmDup.bam.bai,percell_dup_metrics.tsv}'

    input:
    tuple val(meta), path(bams), path(bais), val(cell_ids)

    output:
    tuple val(meta), path("rmdup/*.bam"), path("rmdup/*.bam.bai"), val(cell_ids), emit: bams
    tuple val(meta), path("${meta.id}_rmDup.bam"), path("${meta.id}_rmDup.bam.bai"), emit: merged
    tuple val(meta), path("${meta.id}_percell_dup_metrics.tsv"), emit: metrics

    script:
    def bam_list  = bams instanceof List ? bams : [bams]
    def cid_list  = cell_ids instanceof List ? cell_ids : [cell_ids]
    def bam_lines = bam_list.collect { it.name }.join('\n')
    def cid_lines = cid_list.join('\n')
    def heap_gb   = Math.max(1, ((task.memory.toGiga() * 0.8) / task.cpus).intValue())
    """
    #!/bin/bash
    set -euo pipefail
    mkdir -p rmdup

    cat <<'BAMEOF' > _bam_list.txt
${bam_lines}
BAMEOF
    cat <<'CIDEOF' > _cid_list.txt
${cid_lines}
CIDEOF
    if [ "\$(wc -l < _bam_list.txt)" -ne "\$(wc -l < _cid_list.txt)" ]; then
        echo "ERROR: ${meta.id}: number of BAMs and cell_ids differ" >&2
        exit 1
    fi

    # One JVM per cell, ${task.cpus} in parallel, heap split across them.
    # MarkDuplicates (GATK 4.6) needs a read group on every read; per-cell BAMs often have none (STAR without
    # --outSAMattrRGline), so each cell gets its own (ID = LB = cell index: one library per cell, SM = sample).
    awk '{printf "%06d\\t%s\\n", NR, \$1}' _bam_list.txt > _jobs.tsv
    xargs -P ${task.cpus} -L 1 bash -c '
        samtools addreplacerg -w -r "ID:\$0" -r "SM:${meta.id}" -r "LB:\$0" -r "PL:illumina" \\
            -o rmdup/\$0.rg.bam "\$1" \\
        && gatk --java-options "-Xmx${heap_gb}G" MarkDuplicates \\
            -I rmdup/\$0.rg.bam \\
            -O rmdup/\$0.bam \\
            -M rmdup/\$0.metrics.txt \\
            --REMOVE_DUPLICATES true \\
            --QUIET true \\
            --VERBOSITY WARNING \\
        && rm -f rmdup/\$0.rg.bam \\
        && samtools index rmdup/\$0.bam
    ' < _jobs.tsv

    # Per-cell duplication table (MarkDuplicates metrics row, prefixed by cell_id)
    paste <(cut -f1 _jobs.tsv) _cid_list.txt > _idx_cid.tsv
    FIRST=1
    while IFS=\$'\\t' read -r IDX CID; do
        M=rmdup/\${IDX}.metrics.txt
        if [ "\$FIRST" -eq 1 ]; then
            printf 'sample\\tcell_id\\t%s\\n' "\$(awk '/^LIBRARY/{print; exit}' \$M)" > ${meta.id}_percell_dup_metrics.tsv
            FIRST=0
        fi
        printf '%s\\t%s\\t%s\\n' "${meta.id}" "\$CID" "\$(awk '/^LIBRARY/{getline; print; exit}' \$M)" >> ${meta.id}_percell_dup_metrics.tsv
    done < _idx_cid.tsv
    rm -f rmdup/*.metrics.txt

    if [ "\$(ls rmdup/*.bam | wc -l)" -eq 1 ]; then
        cp rmdup/*.bam ${meta.id}_rmDup.bam
    else
        samtools merge -@ ${task.cpus} ${meta.id}_rmDup.bam rmdup/*.bam
    fi
    samtools index ${meta.id}_rmDup.bam
    """
}


/*
 * UMI dedup per chromosome and GEM well (dedup_per_well_script). Dedup is
 * position-local, so per-chromosome runs equal one pass, in bounded memory.
 * `_other` = every mapped contig not in --main_chroms.
 */
process UMI_DEDUP_CHROM {
    tag { "${meta.id}:${chrom}" }
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), val(chrom)
    path main_chroms

    output:
    tuple val(meta), val(chrom), path("${meta.id}__${chrom}_umiDup.bam"), emit: bam

    script:
    def cell_tag = params.dedup_cell_tag ?: 'CR'
    def out = "${meta.id}__${chrom}_umiDup.bam"
    """
    #!/bin/bash
    set -euo pipefail

    HAS_UB=\$(nf_head.sh 1000 samtools view ${bam} | grep -c "UB:Z:" || true)
    if [ "\$HAS_UB" -eq 0 ]; then
        echo "ERROR: UMI-aware dedup requires UB (corrected UMI) tag but none found." >&2
        echo "Use --dedup coordinate or --dedup none instead." >&2
        exit 1
    fi

    dedup_unit() {
        umi_tools dedup --per-cell --cell-tag=${cell_tag} --umi-tag=UB \\
            --extract-umi-method=tag --method=unique \\
            -I "\$1" -S "\$2" --log="\${2%.bam}.log"
    }
    ${dedup_per_well_script(bam, chrom, main_chroms, cell_tag, task.cpus, out)}
    cat _dd/*.log > ${meta.id}__${chrom}_umi_dedup.log 2>/dev/null || true
    rm -rf _dd
    """
}

// Checks that no read is lost in the merge.
process UMI_DEDUP_GATHER {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/02_dedup" }, mode: 'copy'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}_umiDup.bam"), path("${meta.id}_umiDup.bam.bai"), emit: bam

    script:
    def bl = (bams instanceof List ? bams : [bams]).collect { it.name }.sort().join(' ')
    """
    #!/bin/bash
    set -euo pipefail
    SUM=0
    for b in ${bl}; do SUM=\$((SUM + \$(samtools view -c \$b))); done
    samtools merge -@ ${task.cpus} -c -p -f -o ${meta.id}_umiDup.bam ${bl}
    samtools index -@ ${task.cpus} ${meta.id}_umiDup.bam
    N=\$(samtools view -c ${meta.id}_umiDup.bam)
    [ "\$N" -eq "\$SUM" ] || { echo "ERROR: merged \$N reads != sum of pieces \$SUM" >&2; exit 1; }
    echo "UMI dedup ${meta.id}: \$N reads from \$(echo ${bl} | wc -w) chromosome pieces"
    """
}


/*
 * --dedup none: symlink the input BAM.
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
