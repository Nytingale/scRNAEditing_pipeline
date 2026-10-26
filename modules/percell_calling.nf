/*
 * Per-cell variant calling and AEI
 */

process PERCELL_MPILEUP {
    tag { "${meta.id}" }
    label 'process_high'
    publishDir { "${params.outdir}/${meta.id}/05_percell_variants" }, mode: 'copy'

    input:
    tuple val(meta), path(percell_bams)
    path fasta

    output:
    tuple val(meta), path("percell_variants/*.txt.gz"), emit: variants

    script:
    def qual = params.qual
    """
    mkdir percell_variants

    for BAM in ${percell_bams}; do
        CELL=\$(basename \$BAM .bam)

        bcftools mpileup \\
            -d 100000000 \\
            -q ${qual} -Q ${qual} \\
            -f ${fasta} \\
            \$BAM 2>/dev/null | \\
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP]\\t[%AD]\\n' | \\
        gzip > percell_variants/\${CELL}_variants.txt.gz || true
    done
    """
}


process PERCELL_AEI {
    tag { "${meta.id}" }
    label 'process_high'
    publishDir { "${params.outdir}/${meta.id}/05_percell_aei" }, mode: 'copy'

    input:
    tuple val(meta), path(percell_bams)
    path repeat_annots
    path fasta

    output:
    tuple val(meta), path("${meta.id}_percell_aei.tsv"), emit: aei

    script:
    """
    echo -e "cell\\tA2G_EI\\ttotal_EI\\tpercent_A2G" > ${meta.id}_percell_aei.tsv

    for BAM in ${percell_bams}; do
        CELL=\$(basename \$BAM .bam)
        NREADS=\$(samtools view -c \$BAM)

        if [ "\$NREADS" -lt 300 ]; then
            echo -e "\${CELL}\\tNA\\tNA\\tNA" >> ${meta.id}_percell_aei.tsv
            continue
        fi

        samtools mpileup -l ${repeat_annots} -d 100000 -q 30 -Q 30 -f ${fasta} \$BAM 2>/dev/null | \\
        awk -F'\\t' '
        {
            ref = toupper(\$3)
            if (ref != "A" && ref != "T") next

            bases = \$5
            total_ref = 0
            editing = 0

            n = split(bases, arr, "")
            for (i = 1; i <= n; i++) {
                c = arr[i]
                if (ref == "A") {
                    if (c == "." || c == ",") total_ref++
                    else if (c == "g" || c == "G") editing++
                } else if (ref == "T") {
                    if (c == "." || c == ",") total_ref++
                    else if (c == "c" || c == "C") editing++
                }
            }

            sum_ref += total_ref + editing
            sum_edit += editing
        }
        END {
            if (sum_ref > 0) {
                printf "%s\\t%d\\t%d\\t%.6f\\n", "'\$CELL'", sum_edit, sum_ref, sum_edit/sum_ref*100
            } else {
                printf "%s\\t0\\t0\\tNA\\n", "'\$CELL'"
            }
        }' >> ${meta.id}_percell_aei.tsv
    done
    """
}


/*
 * Per-cell discovery on the pre-split per-cell BAMs (SINTO_SPLIT_PERCELL) with any of
 * the five callers; bin/downstream_filter_percell.R applies the editing filters per
 * cell. Per-cell AEI is PERCELL_AEI (step 6).
 */
process PERCELL_CALLING {
    tag { "${meta.id}:${params.percell_caller}" }
    label 'process_high'
    publishDir { "${params.outdir}/${meta.id}/06_percell_discovery" }, mode: 'copy'

    input:
    tuple val(meta), path(percell_bams)
    path fasta
    path fasta_fai
    path fasta_dict
    path snp_file
    path strand_annots
    path main_chroms

    output:
    tuple val(meta), path("${meta.id}_percell_discovery.tsv"), emit: discovery

    script:
    def caller = params.percell_caller
    def qual   = params.qual
    def gatk_xmx = (task.memory.toGiga() * 0.5).intValue()
    """
    #!/bin/bash
    set -eo pipefail

    SAMPLE="${meta.id}"
    CALLER="${caller}"
    QUAL=${qual}
    FASTA="${fasta}"

    TMPWORK="\${TMPDIR:-/tmp}/percell_\${SAMPLE}_\$\$"
    mkdir -p "\$TMPWORK"

    # Resolve JACUSA2 jar once (jacusa2 branch only).
    JACUSA2_JAR=\$(ls \${CONDA_PREFIX:-}/share/jacusa2-*/JACUSA_v*.jar 2>/dev/null | head -n1)
    if [ -z "\$JACUSA2_JAR" ]; then
        echo "ERROR: no JACUSA_v*.jar under \${CONDA_PREFIX:-<no conda env>}/share/jacusa2-*/ (install jacusa2 into the env)" >&2
        exit 1
    fi

    RAW_TSV="\$TMPWORK/raw_percell_variants.tsv"
    echo -e "cell_id\\tchr\\tpos\\tref\\talt\\tdp\\tad\\ti16" > "\$RAW_TSV"

    TOTAL_CELLS=\$(ls -1 ${percell_bams} 2>/dev/null | wc -l)
    echo "=== \$(date) === Per-cell discovery (\$CALLER) for \$SAMPLE: \$TOTAL_CELLS cells ==="
    CELL_NUM=0

    # Iterate the pre-split, sorted, indexed per-cell BAMs directly.
    for CELL_BAM in ${percell_bams}; do
        BC=\$(basename "\$CELL_BAM" .bam)
        CELL_NUM=\$((CELL_NUM + 1))
        if (( CELL_NUM % 500 == 0 )); then
            echo "  \$(date +%H:%M:%S)  Processed \$CELL_NUM / \$TOTAL_CELLS cells"
        fi
        # Ensure an index exists (SINTO emits sorted BAMs; .bai may not stage).
        [ -f "\${CELL_BAM}.bai" ] || samtools index "\$CELL_BAM" 2>/dev/null || true

        if [ "\$CALLER" = "bcftools" ]; then
            bcftools mpileup -d 100000000 -q \$QUAL -Q \$QUAL -a FORMAT/AD,FORMAT/DP -I -f "\$FASTA" "\$CELL_BAM" 2>/dev/null | \\
            bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP]\\t[%AD]\\t%INFO/I16\\n' | \\
            awk -v cell="\$BC" -v OFS='\\t' '{print cell, \$0}' >> "\$RAW_TSV" || true

        elif [ "\$CALLER" = "reditools3" ]; then
            python -m reditools analyze -r "\$FASTA" -s 0 -q \$QUAL -bq \$QUAL -V \\
                -o "\$TMPWORK/redi.txt" "\$CELL_BAM" 2>/dev/null || true
            if [ -f "\$TMPWORK/redi.txt" ]; then
                awk -v cell="\$BC" -v OFS='\\t' 'NR>1{print cell, \$1, \$2, \$3, \$8, \$5, \$7, "."}' \\
                    "\$TMPWORK/redi.txt" >> "\$RAW_TSV" || true
                rm -f "\$TMPWORK/redi.txt"
            fi

        elif [ "\$CALLER" = "gatk" ]; then
            CELL_RG="\$TMPWORK/rg.bam"
            gatk --java-options "-Xmx${gatk_xmx}G" AddOrReplaceReadGroups \\
                -I "\$CELL_BAM" -O "\$CELL_RG" \\
                -RGID "\$SAMPLE" -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM "\$SAMPLE" 2>/dev/null || true
            samtools index "\$CELL_RG" 2>/dev/null || true
            gatk --java-options "-Xmx${gatk_xmx}G" HaplotypeCaller \\
                -R "\$FASTA" -I "\$CELL_RG" -O "\$TMPWORK/c.vcf" \\
                --dont-use-soft-clipped-bases true 2>/dev/null || true
            if [ -f "\$TMPWORK/c.vcf" ]; then
                grep -v '^#' "\$TMPWORK/c.vcf" | \\
                awk -v cell="\$BC" -v OFS='\\t' '{print cell, \$1, \$2, \$4, \$5, ".", ".", "."}' >> "\$RAW_TSV" || true
                rm -f "\$TMPWORK/c.vcf" "\$TMPWORK/c.vcf.idx"
            fi
            rm -f "\$CELL_RG" "\$CELL_RG.bai"

        elif [ "\$CALLER" = "scallele" ]; then
            ( eval "\$(conda shell.bash hook 2>/dev/null)"; mamba activate scAllele; \\
              scAllele -b "\$CELL_BAM" -g "\$FASTA" -o "\$TMPWORK/sa" \\
                  --run_mode Variant_Caller \\
                  --min-base_quality \$QUAL --min-map_quality \$QUAL 2>/dev/null ) || true
            if [ -f "\$TMPWORK/sa.vcf" ]; then
                scallele_io.py vcf "\$TMPWORK/sa.vcf" 1 2>/dev/null | \\
                awk -v cell="\$BC" -v OFS='\\t' '{print cell, \$1, \$2, \$3, \$4, \$7, \$6, "."}' >> "\$RAW_TSV" || true
            fi
            rm -f "\$TMPWORK/sa"*

        elif [ "\$CALLER" = "jacusa2" ]; then
            java -Xmx4G -jar "\$JACUSA2_JAR" call-1 -R "\$FASTA" -c 1 -m \$QUAL -q \$QUAL \\
                -p 1 -P UNSTRANDED -r "\$TMPWORK/j2.out" "\$CELL_BAM" 2>/dev/null || true
            if [ -f "\$TMPWORK/j2.out" ]; then
                awk -F'\\t' -v cell="\$BC" -v OFS='\\t' '
                    /^##/ {next}
                    /^#/  { for(i=1;i<=NF;i++){h=\$i; sub(/^#/,"",h); col[h]=i} next }
                    {
                        b=\$(col["bases11"]); ref=toupper(\$(col["ref"]));
                        split(b,ct,","); nA=ct[1]; nC=ct[2]; nG=ct[3]; nT=ct[4];
                        dp=nA+nC+nG+nT; if(dp==0) next;
                        split("A C G T",bn," "); cnt[bn[1]]=nA;cnt[bn[2]]=nC;cnt[bn[3]]=nG;cnt[bn[4]]=nT;
                        alt=""; amax=-1;
                        for(k=1;k<=4;k++){ if(bn[k]!=ref && cnt[bn[k]]>amax){amax=cnt[bn[k]]; alt=bn[k]} }
                        if(alt=="" || amax<=0) next;
                        print cell, \$(col["contig"]), \$(col["end"]), ref, alt, dp, amax, "."
                    }' "\$TMPWORK/j2.out" >> "\$RAW_TSV" || true
                rm -f "\$TMPWORK/j2.out"*
            fi
        else
            echo "ERROR: Unknown caller \$CALLER" >&2
            exit 1
        fi
    done

    echo "  Total variant lines: \$(( \$(wc -l < "\$RAW_TSV") - 1 ))"

    # ── Filtering: downstream_filter_percell.R (all cells at once) ────────
    echo "=== \$(date) === Per-cell filtering ==="
    gzip "\$RAW_TSV"

    Rscript ${projectDir}/bin/downstream_filter_percell.R \\
        --input "\${RAW_TSV}.gz" \\
        --caller "\$CALLER" \\
        --sample "\$SAMPLE" \\
        --genome_fasta "\$FASTA" \\
        --snp_rdata "${snp_file}" \\
        --strand_annots "${strand_annots}" \\
        --main_chroms "${main_chroms}" \\
        --min_cov ${params.min_cov} \\
        --min_alt ${params.min_alt ?: 3} \\
        --min_edit_level ${params.min_edit_level ?: 0.01} \\
        --output "${meta.id}_percell_discovery.tsv"

    rm -rf "\$TMPWORK"
    echo "=== \$(date) === Per-cell discovery complete ==="
    """
}
