/*
 * BAM merging modules
 *
 * MERGE_BAMS: Standard merge of all BAMs sharing a sample_id.
 *   Used for sample-level BAM types (sample, sampleRmDup).
 *
 * SELECTIVE_MERGE: Merge only a subset of per-cell BAMs matching an
 *   annotation/filter file. Used for goodCells, annotated, and celltype
 *   BAM types on per-cell platforms (e.g. Smart-seq2) where sinto cannot
 *   work (no CB tags in BAM). The annotation file maps cell_id → group;
 *   only BAMs whose cell_id appears in column 1 of the annotation are merged.
 *   For celltype, produces one BAM per unique group (column 2).
 */

process MERGE_BAMS {
    tag { "${meta.id}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/01_merged" }, mode: 'copy'

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.id}_merged.bam"), path("${meta.id}_merged.bam.bai"), emit: bam

    script:
    def input_bams = bams instanceof List ? bams.join(' ') : bams
    """
    if [ \$(echo "${input_bams}" | wc -w) -eq 1 ]; then
        cp -L ${input_bams} ${meta.id}_merged.bam
    else
        samtools merge -@ ${task.cpus} ${meta.id}_merged.bam ${input_bams}
    fi
    samtools index ${meta.id}_merged.bam
    """
}


/*
 * SELECTIVE_MERGE: merge a subset of per-cell BAMs matching an annotation file.
 *
 * Inputs:
 *   meta         — sample metadata map (must have .id)
 *   bams         — list of per-cell BAMs (staged by Nextflow)
 *   bais         — list of corresponding BAI files
 *   cell_ids     — list of cell identifiers, parallel to bams (same order)
 *   annot_file   — TSV with cell_id in col1 and group in col2
 *   filter_label — output label (e.g. 'goodCells', 'annotated')
 *   mode         — 'filter' (merge all matched into one BAM) or
 *                  'celltype' (one BAM per unique group in col2)
 *
 * Outputs (filter mode):
 *   bam: [meta, <label>.bam, <label>.bam.bai]
 *
 * Outputs (celltype mode):
 *   bam: [meta, <celltype>.bam, <celltype>.bam.bai]  (one per cell type)
 */
process SELECTIVE_MERGE {
    tag { "${meta.id}_${filter_label}" }
    label 'process_medium'
    publishDir { "${params.outdir}/${meta.id}/02_bam_types/${filter_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bams), path(bais), val(cell_ids), path(annot_file), val(filter_label), val(mode)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam

    script:
    def bam_list  = bams instanceof List ? bams : [bams]
    def cid_list  = cell_ids instanceof List ? cell_ids : [cell_ids]
    def bam_lines = bam_list.collect { it.name }.join('\n')
    def cid_lines = cid_list.join('\n')
    """
    # Build cell_id → BAM filename mapping
    cat <<'BAMEOF' > _bam_list.txt
${bam_lines}
BAMEOF
    cat <<'CIDEOF' > _cid_list.txt
${cid_lines}
CIDEOF
    paste _cid_list.txt _bam_list.txt > _cellid_to_bam.tsv
    NBAMS=\$(wc -l < _cellid_to_bam.tsv)

    if [ "${mode}" = "filter" ]; then
        # Filter mode: merge all BAMs whose cell_id appears in col1 of annot
        cut -f1 ${annot_file} | sort -u > _wanted_cells.txt
        BAMS_TO_MERGE=""
        while IFS=\$'\\t' read -r CID BAM; do
            if grep -qFx "\$CID" _wanted_cells.txt; then
                BAMS_TO_MERGE="\$BAMS_TO_MERGE \$BAM"
            fi
        done < _cellid_to_bam.tsv

        N_MATCH=\$(echo \$BAMS_TO_MERGE | wc -w)
        if [ "\$N_MATCH" -eq 0 ]; then
            echo "ERROR: No cell_ids from annotation matched any BAM cell_ids" >&2
            echo "  Annotation file: ${annot_file}" >&2
            echo "  First 5 annotation cell_ids:" >&2
            head -5 _wanted_cells.txt >&2
            echo "  First 5 BAM cell_ids:" >&2
            head -5 _cid_list.txt >&2
            exit 1
        fi
        echo "Merging \$N_MATCH / \$NBAMS BAMs for ${filter_label}"

        if [ "\$N_MATCH" -eq 1 ]; then
            cp -L \$BAMS_TO_MERGE ${meta.id}_${filter_label}.bam
        else
            samtools merge -@ ${task.cpus} ${meta.id}_${filter_label}.bam \$BAMS_TO_MERGE
        fi
        samtools index ${meta.id}_${filter_label}.bam

    elif [ "${mode}" = "celltype" ]; then
        # Celltype mode: one merged BAM per unique group (col2 of annot)
        # Clean annotation: sanitize spaces in cell type names
        awk -F'\\t' '{
            gsub(/[[:space:]]+\$/, "", \$1);
            gsub(/[[:space:]]+\$/, "", \$2);
            gsub(/ /, "_", \$2);
            if (\$1 != "" && \$2 != "") print \$1"\\t"\$2
        }' ${annot_file} > _clean_annot.tsv

        cut -f2 _clean_annot.tsv | sort -u > _celltypes.txt
        while read -r CT; do
            # Get cell_ids belonging to this cell type
            awk -F'\\t' -v ct="\$CT" '\$2 == ct {print \$1}' _clean_annot.tsv | sort -u > _ct_cells.txt
            BAMS_TO_MERGE=""
            while IFS=\$'\\t' read -r CID BAM; do
                if grep -qFx "\$CID" _ct_cells.txt; then
                    BAMS_TO_MERGE="\$BAMS_TO_MERGE \$BAM"
                fi
            done < _cellid_to_bam.tsv

            N_CT=\$(echo \$BAMS_TO_MERGE | wc -w)
            if [ "\$N_CT" -gt 0 ]; then
                echo "Celltype \$CT: merging \$N_CT BAMs"
                if [ "\$N_CT" -eq 1 ]; then
                    cp -L \$BAMS_TO_MERGE \${CT}.bam
                else
                    samtools merge -@ ${task.cpus} \${CT}.bam \$BAMS_TO_MERGE
                fi
                samtools index \${CT}.bam
            else
                echo "WARNING: No BAMs found for celltype \$CT (skipping)"
            fi
        done < _celltypes.txt
    fi
    """
}
