/*
 * sinto barcode filtering and cell-type splitting. SINTO_EXTRACT_QC is a generic
 * barcode filter, aliased in main.nf as SINTO_GOOD_CELLS* and SINTO_ANNOTATED*.
 */

process SINTO_EXTRACT_QC {
    tag { "${meta.id}:${filter_label}" }
    label 'process_medium'
    // Time: nextflow.config. 'finish' keeps a failure here from killing caller
    // tasks already in flight.
    errorStrategy { task.exitStatus in [137, 140, 143, 104, 71] && task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries    2
    publishDir { "${params.outdir}/${meta.id}/02_barcode_filtered/${filter_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(barcode_source), val(filter_label)

    output:
    tuple val(meta), path("${meta.id}_${filter_label}.bam"), path("${meta.id}_${filter_label}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_${filter_label}_barcode_coverage.txt"), emit: coverage

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Extract unique barcodes from the source file (col 1), create 2-column format for sinto
    # Sanitize: strip trailing whitespace/carriage returns from barcodes
    awk -F'\\t' '{gsub(/[[:space:]]+\$/, "", \$1); if (\$1 != "") print \$1"\\tfiltered"}' ${barcode_source} | sort -u -k1,1 > _filter_barcodes.txt
    N_FILTER_BC=\$(wc -l < _filter_barcodes.txt)

    if [ "\$N_FILTER_BC" -eq 0 ]; then
        echo "ERROR: No barcodes extracted from ${barcode_source}" >&2
        exit 1
    fi

    # Spot-check: do any filter barcodes match CB tags in the BAM?
    nf_head.sh 100000 samtools view ${bam} | grep -oP 'CB:Z:\\K[^\\t]+' | sort -u > _bam_sample_bcs.txt || true
    N_OVERLAP=\$(comm -12 <(cut -f1 _filter_barcodes.txt | sort) _bam_sample_bcs.txt | wc -l)

    if [ "\$N_OVERLAP" -eq 0 ]; then
        echo "ERROR: Zero barcode overlap between filter list and BAM CB tags." >&2
        echo "  Filter barcodes (first 5):" >&2
        head -5 _filter_barcodes.txt | cut -f1 >&2
        echo "  BAM barcodes (first 5):" >&2
        head -5 _bam_sample_bcs.txt >&2
        echo "  Common cause: barcode suffix mismatch (e.g. '-1' in CellRanger vs no suffix in STARsolo)." >&2
        echo "  Check that barcode formats match between your annotation file and BAM CB tags." >&2
        exit 1
    fi
    echo "Barcode spot-check: \$N_OVERLAP / \$N_FILTER_BC filter barcodes found in BAM sample"
    rm -f _bam_sample_bcs.txt

    # Estimate total unique barcodes in the BAM (CB tag) — from the same sample
    N_BAM_BC=\$(nf_head.sh 1000000 samtools view ${bam} | grep -oP 'CB:Z:\\K[^\\t]+' | sort -u | grep -cv '^-\$' || echo 0)

    sinto filterbarcodes \\
        -b ${bam} \\
        -c _filter_barcodes.txt \\
        -p ${task.cpus} \\
        --outdir .

    mv filtered.bam ${meta.id}_${filter_label}.bam 2>/dev/null || \\
        samtools view -h -b ${bam} -D CB:<(cut -f1 _filter_barcodes.txt) > ${meta.id}_${filter_label}.bam

    samtools index ${meta.id}_${filter_label}.bam
    N_RETAINED_READS=\$(samtools view -c ${meta.id}_${filter_label}.bam)
    N_TOTAL_READS=\$(samtools view -c ${bam})

    if [ "\$N_RETAINED_READS" -eq 0 ]; then
        echo "ERROR: Filtered BAM has 0 reads despite \$N_FILTER_BC barcodes in filter list." >&2
        echo "  This likely means barcode formats don't match (suffix mismatch)." >&2
        exit 1
    fi

    echo -e "sample\\tfilter_label\\tfilter_barcodes\\tbam_barcodes\\tbarcode_coverage\\tfiltered_reads\\ttotal_reads\\tread_retention" > ${meta.id}_${filter_label}_barcode_coverage.txt
    echo -e "${meta.id}\\t${filter_label}\\t\${N_FILTER_BC}\\t\${N_BAM_BC}\\t\$(echo "scale=4; \$N_FILTER_BC / (\$N_BAM_BC + 0.001)" | bc)\\t\${N_RETAINED_READS}\\t\${N_TOTAL_READS}\\t\$(echo "scale=4; \$N_RETAINED_READS / (\$N_TOTAL_READS + 0.001)" | bc)" >> ${meta.id}_${filter_label}_barcode_coverage.txt

    rm -f _filter_barcodes.txt
    """
}


process SINTO_EXTRACT_CELLTYPE {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/03_celltype_bams" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(celltype_annot)

    output:
    tuple val(meta), path("*.celltype.bam"), path("*.celltype.bam.bai"), emit: celltype_bams

    script:
    // celltype_annot is a 2-column TSV: barcode \t celltype
    def min_ct = (params.min_celltype_cells ?: 0) as int
    """
    #!/bin/bash
    set -euo pipefail

    # Sanitize celltype annotation: replace spaces in col 2 with underscores,
    # strip trailing whitespace/carriage returns from both columns
    awk -F'\\t' '{
        gsub(/[[:space:]]+\$/, "", \$1);
        gsub(/[[:space:]]+\$/, "", \$2);
        gsub(/ /, "_", \$2);
        if (\$1 != "" && \$2 != "") print \$1"\\t"\$2
    }' ${celltype_annot} > _clean_annot.tsv

    # Optional min-cell filter: keep only cell types with >= min_celltype_cells
    # cells (barcodes). Used by the celltype-as-sample discovery mode so tiny,
    # under-powered cell types are not run through the full caller matrix.
    MIN_CT=${min_ct}
    if [ "\$MIN_CT" -gt 0 ]; then
        cut -f2 _clean_annot.tsv | sort | uniq -c | awk -v m="\$MIN_CT" '\$1>=m{print \$2}' > _keep_ct.txt
        awk -F'\\t' 'NR==FNR{k[\$0]=1; next} (\$2 in k)' _keep_ct.txt _clean_annot.tsv > _clean_annot.flt.tsv
        mv _clean_annot.flt.tsv _clean_annot.tsv
        echo "Cell-type min-cell filter (>= \$MIN_CT): kept \$(wc -l < _keep_ct.txt) cell types"
    fi
    N_ANNOT_BC=\$(wc -l < _clean_annot.tsv)

    if [ "\$N_ANNOT_BC" -eq 0 ]; then
        echo "ERROR: No valid barcode-celltype pairs in ${celltype_annot}" >&2
        exit 1
    fi

    # GUARD (read misassignment): a single annotation must map each barcode to ONE
    # cell type. A barcode appearing under >1 cell type means the annotation pools
    # multiple samples/lanes with reused barcodes (10x barcodes are NOT unique
    # across samples) — sinto would then assign that barcode's reads to a single
    # arbitrary cell type. Fail loud rather than silently misassign. (The
    # pipeline's split_by_celltype mode avoids this entirely by splitting
    # per-sample BEFORE merge; this guard protects any other/legacy usage.)
    DUP_BC=\$(cut -f1 _clean_annot.tsv | sort | uniq -d | head -1)
    if [ -n "\$DUP_BC" ]; then
        N_DUP=\$(cut -f1 _clean_annot.tsv | sort | uniq -d | wc -l)
        echo "ERROR: \$N_DUP barcode(s) map to >1 cell type in ${celltype_annot} (e.g. '\$DUP_BC')." >&2
        echo "  This annotation pools samples/lanes with REUSED barcodes; a raw-CB split" >&2
        echo "  would MISASSIGN reads. Provide per-sample annotations and split per sample" >&2
        echo "  before merging (pipeline --split_by_celltype does this automatically)." >&2
        exit 1
    fi

    # Spot-check: do annotation barcodes match CB tags in the BAM?
    nf_head.sh 100000 samtools view ${bam} | grep -oP 'CB:Z:\\K[^\\t]+' | sort -u > _bam_sample_bcs.txt || true
    N_OVERLAP=\$(comm -12 <(cut -f1 _clean_annot.tsv | sort -u) _bam_sample_bcs.txt | wc -l)

    if [ "\$N_OVERLAP" -eq 0 ]; then
        echo "ERROR: Zero barcode overlap between celltype annotation and BAM CB tags." >&2
        echo "  Annotation barcodes (first 5):" >&2
        head -5 _clean_annot.tsv | cut -f1 >&2
        echo "  BAM barcodes (first 5):" >&2
        head -5 _bam_sample_bcs.txt >&2
        echo "  Common cause: barcode suffix mismatch (e.g. '-1' in CellRanger vs no suffix in STARsolo)." >&2
        echo "  Check that barcode formats match between your annotation file and BAM CB tags." >&2
        exit 1
    fi
    echo "Barcode spot-check: \$N_OVERLAP / \$N_ANNOT_BC annotation barcodes found in BAM sample"
    rm -f _bam_sample_bcs.txt

    sinto filterbarcodes \\
        -b ${bam} \\
        -c _clean_annot.tsv \\
        --outdir . \\
        -p ${task.cpus}

    # sinto names output as <celltype>.bam; rename for clarity and index
    # Skip the input BAM and annotation file (staged by Nextflow in the same directory)
    INPUT_BAM_NAME=\$(basename ${bam})
    N_CT_BAMS=0
    for f in *.bam; do
        [ "\$f" = "\$INPUT_BAM_NAME" ] && continue
        ct=\$(basename \$f .bam)
        mv \$f ${meta.id}_\${ct}.celltype.bam
        samtools index ${meta.id}_\${ct}.celltype.bam
        N_CT_BAMS=\$((N_CT_BAMS + 1))
    done

    if [ "\$N_CT_BAMS" -eq 0 ]; then
        echo "ERROR: sinto produced no celltype BAMs despite \$N_ANNOT_BC annotated barcodes." >&2
        echo "  This likely means barcode formats don't match (suffix mismatch)." >&2
        exit 1
    fi
    echo "Produced \$N_CT_BAMS celltype BAMs from \$N_ANNOT_BC annotated barcodes"

    rm -f _clean_annot.tsv
    """
}


process SINTO_SPLIT_PERCELL {
    tag { "${meta.id}" }
    label 'process_high'
    // cpus / memory: nextflow.config (sinto's per-cell file handles scale RAM with cpus).
    errorStrategy { task.exitStatus in [137, 140, 143, 104, 71] ? 'retry' : 'terminate' }
    maxRetries 2
    publishDir { "${params.outdir}/${meta.id}/04_percell_bams" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes)

    output:
    tuple val(meta), path("percell/*.bam"), emit: percell_bams

    script:
    // Restrict to the called barcodes when given: the sample BAM holds many
    // background barcodes.
    """
    #!/bin/bash
    set -euo pipefail
    mkdir percell

    if [ "${barcodes}" != "NO_FILE" ] && [ -s "${barcodes}" ]; then
        cut -f1 "${barcodes}" | sed 's/[[:space:]]*\$//' | grep -v '^\$' | sort -u \\
            | awk '{print \$1"\\t"\$1}' > _percell_barcodes.txt
    else
        samtools view ${bam} | grep -oP 'CB:Z:\\K[^\\t]+' | sort -u | awk '{print \$1"\\t"\$1}' > _percell_barcodes.txt
    fi

    sinto filterbarcodes \\
        -b ${bam} \\
        -c _percell_barcodes.txt \\
        --outdir percell \\
        -p ${task.cpus}

    rm -f _percell_barcodes.txt

    # Sort and index each per-cell BAM
    for f in percell/*.bam; do
        samtools sort -o \${f}.sorted \$f
        mv \${f}.sorted \$f
        samtools index \$f
    done
    """
}
