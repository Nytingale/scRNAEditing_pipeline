#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * scRNAEditing: Benchmarking pipeline for A-to-I RNA editing in scRNA-seq
 *
 * Ng & Tan, 2026
 * https://github.com/Nytingale/scRNAEditing
 */

params.input           = null    // Samplesheet CSV
params.outdir          = 'results'
params.genome_fasta    = null
params.genome_gtf      = null
params.star_index      = null
params.snp_file        = null    // Known SNPs: RData (GRanges, preferred) or BED
params.known_sites_vcf = null    // GATK BQSR known sites (null = no recalibration)
params.strand_annots   = null
params.alu_annots      = null
params.rediportal      = null    // REDIportal VCF; when set, writes summary/rediportal_overlap.tsv
params.repeat_annots   = null    // Repeat annotations for AEI
params.refseq_bed      = null    // RefSeq BED for strandedness inference
params.main_chroms     = null    // One chromosome per line (null = all from BAM header)
params.remap           = false   // Re-align with STARsolo
params.tso_clip        = true    // --remap: 10x 3' -> STAR --clipAdapterType CellRanger4 (5' TSO + 3' poly(A)); smartseq2 -> cutadapt 5' TSO trim on both mates
params.tso_seq         = 'AAGCAGTGGTATCAACGCAGAGTACATGGG'
params.tso_min_overlap = 5       // shortest TSO 3' end (ATGGG) trimmed at a read start
params.tso_min_length  = 20      // drop pairs with a mate shorter than this after trimming
params.cutadapt_conda  = "${projectDir}/envs/cutadapt.yml"  // or an existing env prefix
params.dedup           = 'none,coordinate'  // Any of: none, coordinate, umi
params.dedup_cell_tag  = 'CR'          // 'CR' (raw, recommended) or 'CB' (corrected)
params.callers         = 'bcftools,reditools3,gatk,jacusa2'  // also: scallele, cellsnp
params.bam_types       = 'sample,sampleRmDup,sampleUmiDup,recal,recalq20,goodCells,goodCellsRmDup,goodCellsUmiDup,annotated,annotatedRmDup,annotatedUmiDup,celltype,celltypeRmDup,celltypeUmiDup'
// Each cell type becomes its own sample (id = <sample>__<celltype>) and runs the
// full sample-level matrix. Use sample-level bam_types only.
params.split_by_celltype   = false
params.min_celltype_cells  = 0   // split_by_celltype: drop cell types with fewer cells (0 = keep all)
params.qual            = 30      // Base/mapping quality threshold
params.bcftools_baq    = true    // true = bcftools default (partial BAQ); false = -B
params.adaptive_windows = '0,20,40,60,80,100,120,140,160,180,200'
// Remove sites with another mismatch type within w1 nt (w1 searched over adaptive_windows).
params.mismatch_window_filter = false
// Keep only sites with another site within w2 nt (off since v8).
params.clustered_filter = false
params.min_cov         = 10
params.strand_modes    = 'annotation'  // Any of: annotation, hybrid, caller
params.strand_ambiguity_threshold = 0.3  // hybrid mode: read-strand ambiguity cutoff

params.percell_analyses      = true    // Step 6; false avoids the memory-heavy SINTO_SPLIT_PERCELL
params.percell_calling       = false   // Step 7
params.percell_aei           = true
params.percell_caller        = null    // Required with percell_calling
params.percell_bam_type      = null    // Required with percell_calling
params.percell_chunk_size    = 100     // Cells per RNAEditingIndexer chunk
params.rnaeditingindexer_path = null
params.aei_genome            = null    // e.g. 'hg38'
params.aei_regions_bed       = null
params.aei_snps_bed          = null
params.aei_refseq_bed        = null
params.aei_expression_bed    = null
params.min_alt               = 3       // alt_rule 'count'
params.min_edit_level        = 0.01
// 'count' (cov_alt >= min_alt) or 'binom' (alt count above the sequencing error rate).
// min_alt, not min_edit_level, is the binding sensitivity threshold.
params.alt_rule              = 'count'
params.seq_error_rate        = 0.001   // binom only
params.alt_binom_alpha       = 0.05    // binom only
params.min_alt_floor         = 2       // binom only
params.jacusa2_min_score     = 1.25
params.jacusa2_scatter       = false   // one JACUSA2 task per chromosome
// Software (envs/, containers/). conda_env and container are read in nextflow.config.
params.conda_env             = null    // existing env prefix; else $SCRNA_CONDA_ENV, else envs/scRNAEditing.linux-64.lock.txt
params.container             = null    // default ghcr.io/nytingale/scrnaediting:<version>
params.scallele_conda        = "${projectDir}/envs/scAllele.yml"  // or an existing env prefix
params.scallele_run_mode     = 'Variant_Caller'  // Variant_Caller | Full | Training
params.cellsnp_bin           = 'cellsnp-lite'    // conda: path to the cellsnp env's binary
// 'off' = dbSNP only; 'on' = dbSNP + the caller's per-sample blacklist (bcftools, gatk only).
params.snp_subtract_modes    = 'off,on'
params.remove_100pct_editing = true
params.save_intermediates    = 0       // Save intermediate tables from step N on (0 = off)
// Add pileup-based editing level and coverage columns (comparable across callers).
params.requant_editing       = true
params.requant_min_bq        = 30
params.requant_min_mq        = 30
params.requant_max_depth     = 1000000
params.mito_chr              = null    // null = auto-detect from GTF
params.qc_min_features       = 200     // goodCells QC of unannotated cells
params.qc_max_features       = 6000
params.qc_max_pct_mt         = 20
params.qc_min_counts         = 500

include { VALIDATE_BAMS        } from './modules/bam_validate'
include { VALIDATE_FASTQ_READS } from './modules/fastq_validate'
include { STARSOLO_ALIGN       } from './modules/starsolo'
include { STAR_GENOMEGENERATE  } from './modules/starsolo'
include { CUTADAPT_TSO         } from './modules/tso_trim'
include { ADD_MD_TAG           } from './modules/bam_validate'
include { FILTER_UNIQUELY_MAPPED } from './modules/bam_validate'
include { MERGE_BAMS           } from './modules/bam_merge'
include { MERGE_BAMS as MERGE_CELLTYPE_BAMS } from './modules/bam_merge'
include { SELECTIVE_MERGE      } from './modules/bam_merge'
include { SELECTIVE_MERGE as SELECTIVE_MERGE_ANNOT   } from './modules/bam_merge'
include { SELECTIVE_MERGE as SELECTIVE_MERGE_CT      } from './modules/bam_merge'
include { SELECTIVE_MERGE as SELECTIVE_MERGE_RMDUP          } from './modules/bam_merge'
include { SELECTIVE_MERGE as SELECTIVE_MERGE_ANNOT_RMDUP    } from './modules/bam_merge'
include { SELECTIVE_MERGE as SELECTIVE_MERGE_CT_RMDUP       } from './modules/bam_merge'
include { SINTO_EXTRACT_QC     as SINTO_GOOD_CELLS           } from './modules/celltype_split'
include { SINTO_EXTRACT_QC     as SINTO_GOOD_CELLS_RMDUP     } from './modules/celltype_split'
include { SINTO_EXTRACT_QC     as SINTO_GOOD_CELLS_UMIDUP   } from './modules/celltype_split'
include { SINTO_EXTRACT_QC     as SINTO_ANNOTATED              } from './modules/celltype_split'
include { SINTO_EXTRACT_QC     as SINTO_ANNOTATED_RMDUP        } from './modules/celltype_split'
include { SINTO_EXTRACT_QC     as SINTO_ANNOTATED_UMIDDUP      } from './modules/celltype_split'
include { SINTO_EXTRACT_CELLTYPE } from './modules/celltype_split'
include { SINTO_EXTRACT_CELLTYPE as SINTO_CELLTYPE_PERLANE } from './modules/celltype_split'
include { SINTO_EXTRACT_CELLTYPE as SINTO_CELLTYPE_RMDUP  } from './modules/celltype_split'
include { SINTO_EXTRACT_CELLTYPE as SINTO_CELLTYPE_UMIDDUP } from './modules/celltype_split'
include { SINTO_SPLIT_PERCELL  } from './modules/celltype_split'
include { SINTO_SPLIT_PERCELL as SINTO_SPLIT_PERCELL_DISC } from './modules/celltype_split'
include { COMPUTE_GOOD_CELLS   } from './modules/goodcells'
include { MARKDUP_CHROM        } from './modules/dedup'
include { MARKDUP_GATHER       } from './modules/dedup'
include { UMI_DEDUP_CHROM      } from './modules/dedup'
include { UMI_DEDUP_GATHER     } from './modules/dedup'
include { MARKDUP_PERCELL_SS2  } from './modules/dedup'
include { GATK_SPLITNCIGAR     } from './modules/gatk_preprocess'
include { GATK_BQSR            } from './modules/gatk_preprocess'
include { INFER_STRANDEDNESS   } from './modules/strandedness'
include { BCFTOOLS_MPILEUP     } from './modules/variant_calling'
include { REDITOOLS3           } from './modules/variant_calling'
include { GATK_HAPLOTYPECALLER } from './modules/variant_calling'
include { JACUSA2_CALL1        } from './modules/variant_calling'
include { JACUSA2_CHROM        } from './modules/variant_calling'
include { JACUSA2_GATHER       } from './modules/variant_calling'
include { CELLSNP_LITE         } from './modules/variant_calling'
include { SCALLELE_CHROM       } from './modules/variant_calling'
include { SCALLELE_GATHER      } from './modules/variant_calling'
include { DOWNSTREAM_R         } from './modules/downstream_r'
include { DOWNSTREAM_R_SCALLELE } from './modules/downstream_r'
include { REQUANT_PILEUP       } from './modules/requant_pileup'
include { COLLATE_RESULTS      } from './modules/collate_results'
include { REDIPORTAL_OVERLAP   } from './modules/rediportal_overlap'
include { PERCELL_MPILEUP      } from './modules/percell_calling'
include { PERCELL_AEI          } from './modules/percell_calling'
include { PERCELL_CALLING      } from './modules/percell_calling'

// `--flag false` arrives as the string "false", which is truthy in Groovy.
def toBool(v) {
    if (v instanceof Boolean) return v
    return (v?.toString()?.trim()?.toLowerCase() in ['true', '1', 'yes', 'y', 't', 'on'])
}

def findSoloFiltered(soloDir) {
    if (soloDir == null) return null
    def root = soloDir.toString()
    def cands = [
        root,
        "${root}/GeneFull/filtered",
        "${root}/Gene/filtered",
        "${root}/GeneFull_Ex50pAS/filtered",
        "${root}/filtered_feature_bc_matrix",
        "${root}/outs/filtered_feature_bc_matrix"
    ]
    def hit = cands.find { c ->
        file("${c}/matrix.mtx").exists() || file("${c}/matrix.mtx.gz").exists()
    }
    return hit ? file(hit) : null
}

def findSoloFromBam(bamFile) {
    if (bamFile == null) return null
    def prefix = bamFile.toString().replaceAll('_Aligned.*', '')
    def solo = file("${prefix}_Solo.out")
    return solo.exists() ? findSoloFiltered(solo) : null
}

// Outputs are named and joined by meta.id, so each cell type needs its own id;
// parent_id keeps the sample id for sample-keyed joins (strandedness).
def celltypeMeta(meta, ct) {
    meta + [ id: "${meta.id}__${ct}".toString(), parent_id: meta.id, celltype: ct ]
}

def flattenCelltype(celltype_out, label) {
    celltype_out.flatMap { meta, bam_list, bai_list ->
        def bams = bam_list instanceof List ? bam_list : [bam_list]
        def bais = bai_list instanceof List ? bai_list : [bai_list]
        [bams, bais].transpose().collect { b, bi ->
            def bn = b.baseName.replaceAll('\\.celltype$', '')
            def ct = bn.replaceFirst('^' + java.util.regex.Pattern.quote(meta.id) + '_', '')
            [ celltypeMeta(meta, ct), b, bi, label ]
        }
    }
}

def flattenSelectiveCelltype(selmerge_out, label) {
    selmerge_out.flatMap { meta, bam_files, bai_files ->
        def bams = bam_files instanceof List ? bam_files : [bam_files]
        def bais = bai_files instanceof List ? bai_files : [bai_files]
        [bams, bais].transpose().collect { b, bi ->
            [ celltypeMeta(meta, b.name.replaceAll('\\.bam$', '')), b, bi, label ]
        }
    }
}

// ── Main workflow ────────────────────────────────────────────────────────────
workflow {

    // ── Validate inputs ─────────────────────────────────────────────────
    if (!params.input)       { error "Please provide --input <samplesheet.csv>" }
    if (!params.genome_fasta){ error "Please provide --genome_fasta <ref.fa>" }
    if (!params.genome_gtf)  { error "Please provide --genome_gtf <genes.gtf>" }
    def valid_dedup = ['none', 'coordinate', 'umi'] as Set
    def requested_dedup = params.dedup.tokenize(',')

    if (toBool(params.percell_calling)) {
        if (!params.percell_caller)   { error "--percell_calling requires --percell_caller <bcftools|reditools3|gatk|jacusa2|scallele>" }
        if (!params.percell_bam_type) { error "--percell_calling requires --percell_bam_type <sample|sampleRmDup|annotated|...>" }
    }
    // ── Parse samplesheet (columns: SAMPLESHEET_SPEC.md) ────────────────
    // Optional columns not in the spec: count_matrix (10x matrix dir; default
    // <bam_prefix>_Solo.out/GeneFull/filtered), qc_barcodes (precomputed goodCells
    // list), rmdup_bam / umidup_bam (prebuilt dedup BAMs, step 3b).
    // Rows sharing a sample_id are lanes of ONE sample and are merged before any
    // analysis. Distinct samples need distinct sample_ids: 10x barcodes are reused
    // across samples, so splitting a merged multi-sample BAM by barcode misassigns
    // reads. To pool cell types across samples use --split_by_celltype.

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def lane_id = row.bam ? file(row.bam).baseName.replaceAll('_Aligned.*', '') :
                          row.fastq_r1 ? file(row.fastq_r1).baseName.replaceAll('_[12]\\.fastq.*', '') :
                          row.sample_id
            def meta = [
                id:        row.sample_id,
                lane_id:   lane_id,
                platform:  row.platform ?: 'unknown',
            ]
            def bam             = row.bam             ? file(row.bam, checkIfExists: true) : null
            def celltype_annot  = file(row.celltype_annot, checkIfExists: true)
            def called_barcodes = row.called_barcodes ? file(row.called_barcodes, checkIfExists: true) : null
            def fastq_r1        = row.fastq_r1        ? file(row.fastq_r1) : null
            def fastq_r2        = row.fastq_r2        ? file(row.fastq_r2) : null
            def cell_id         = row.cell_id          ?: null

            [ meta, bam, celltype_annot, called_barcodes, fastq_r1, fastq_r2, cell_id ]
        }
        .set { ch_input }

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def bam = row.bam ? file(row.bam) : null
            def matrix = row.count_matrix ?
                findSoloFiltered(file(row.count_matrix, checkIfExists: true)) :
                findSoloFromBam(bam)
            def annot  = file(row.celltype_annot, checkIfExists: true)
            def called = row.called_barcodes ? file(row.called_barcodes, checkIfExists: true) : null
            def qc_bc  = row.qc_barcodes ? file(row.qc_barcodes, checkIfExists: true) : null
            [ row.sample_id, row.platform ?: 'unknown', annot, called, matrix, qc_bc, row.cell_id ?: null ]
        }
        .groupTuple(by: 0)
        .map { sid, platforms, annots, calleds, matrices, qcs, cell_ids ->
            def meta = [id: sid, platform: platforms[0]]
            def mats = matrices.findAll { it != null }.collect { it.toString() }.unique()
            def cells = cell_ids.findAll { it != null }
            def called_f = calleds.find { it != null }
            def qc_f = qcs.find { it != null }
            [
                meta,
                annots[0],
                called_f ? called_f.toString() : '',
                qc_f ? qc_f.toString() : '',
                mats ? mats.join(',') : '',
                cells
            ]
        }
        .set { ch_qc_from_sheet }
    def invalid_dedup = requested_dedup.findAll { !(it in valid_dedup) }
    if (invalid_dedup) {
        error "--dedup values must be from: none, coordinate, umi (got invalid: ${invalid_dedup.join(', ')})"
    }

    genome_fasta = file(params.genome_fasta, checkIfExists: true)
    genome_gtf   = file(params.genome_gtf, checkIfExists: true)

    // ── Step 0: Prepare STAR index if remapping ─────────────────────────
    if (toBool(params.remap)) {
        if (params.star_index) {
            ch_star_index = Channel.value(file(params.star_index))
        } else {
            STAR_GENOMEGENERATE(genome_fasta, genome_gtf)
            ch_star_index = STAR_GENOMEGENERATE.out.index
        }
    }

    // ── Step 1: Get per-lane BAMs, then merge by sample_id ────────────
    if (toBool(params.remap)) {
        ch_to_remap = ch_input
            .filter { meta, bam, annot, called_bc, r1, r2, cell_id -> r1 != null && r2 != null }
            .map { meta, bam, annot, called_bc, r1, r2, cell_id -> [ meta, annot, called_bc, r1, r2, cell_id ] }

        VALIDATE_FASTQ_READS(ch_to_remap.map { meta, annot, called_bc, r1, r2, cell_id -> [ meta, r1, r2 ] })

        ch_validated = VALIDATE_FASTQ_READS.out.validated
            .map { meta, bc_fq, cdna_fq -> [ meta.lane_id ?: meta.id, meta, bc_fq, cdna_fq ] }
            .combine(
                ch_to_remap.map { meta, annot, called_bc, r1, r2, cell_id -> [ meta.lane_id ?: meta.id, annot ] },
                by: 0
            )
            .map { lid, meta, bc_fq, cdna_fq, annot -> [ meta, annot, bc_fq, cdna_fq ] }
            .branch { meta, annot, r1, r2 ->
                tso_trim: meta.platform == 'smartseq2' && toBool(params.tso_clip)
                other: true
            }

        CUTADAPT_TSO(ch_validated.tso_trim)
        ch_align_in = ch_validated.other.mix(CUTADAPT_TSO.out.reads)

        STARSOLO_ALIGN(ch_align_in, ch_star_index, genome_fasta, genome_gtf)
        ch_lane_bams = STARSOLO_ALIGN.out.bam
            .map { meta, bam, bai -> [ meta.lane_id, meta, bam, bai ] }
            .combine(
                ch_input.map { meta, bam, annot, called_bc, r1, r2, cell_id -> [ meta.lane_id, annot, called_bc, cell_id ] },
                by: 0
            )
            .map { lane_id, meta, bam, bai, annot, called_bc, cell_id -> [ meta, bam, bai, annot, called_bc, cell_id ] }
    } else {
        ch_lane_input = ch_input
            .filter { meta, bam, annot, called_bc, r1, r2, cell_id -> bam != null }
            .map { meta, bam, annot, called_bc, r1, r2, cell_id -> [ meta, bam, annot, called_bc, cell_id ] }

        VALIDATE_BAMS(ch_lane_input.map { meta, bam, annot, called_bc, cell_id -> [ meta, bam ] })

        ADD_MD_TAG(
            ch_lane_input.map { meta, bam, annot, called_bc, cell_id -> [ meta, bam ] },
            genome_fasta
        )

        FILTER_UNIQUELY_MAPPED(ADD_MD_TAG.out.bam)

        ch_lane_bams = FILTER_UNIQUELY_MAPPED.out.bam
            .map { meta, bam, bai -> [ meta.lane_id, meta, bam, bai ] }
            .combine(
                ch_lane_input.map { meta, bam, annot, called_bc, cell_id -> [ meta.lane_id, annot, called_bc, cell_id ] },
                by: 0
            )
            .map { lane_id, meta, bam, bai, annot, called_bc, cell_id -> [ meta, bam, bai, annot, called_bc, cell_id ] }
    }

    // -> [sample_meta, [bams], [bais], annot, called_bc, [cell_ids] or null]
    ch_for_merge = ch_lane_bams
        .map { meta, bam, bai, annot, called_bc, cell_id ->
            [ meta.id, meta, bam, bai, annot, called_bc, cell_id ]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, bams, bais, annots, called_bcs, cell_ids ->
            def sample_meta = [id: sample_id, platform: metas[0].platform]
            def annot = annots[0]
            // Only the first annotation is used on the merged BAM; differing ones
            // usually mean distinct samples were pooled under one sample_id.
            if (!toBool(params.split_by_celltype)) {
                def distinct_annots = annots.findAll { it != null }.collect { it.toString() }.unique()
                if (distinct_annots.size() > 1) {
                    log.warn "sample_id '${sample_id}' has ${distinct_annots.size()} different celltype_annot files across its rows; only the first is used for barcode-based tiers. Rows sharing a sample_id should be LANES of ONE sample (same annotation). For distinct samples use distinct sample_ids; to pool cell types across samples use --split_by_celltype."
                }
            }
            def called_bc = called_bcs.find { it != null }
            def has_cell_ids = cell_ids.any { it != null }
            [ sample_meta, bams, bais, annot, called_bc, has_cell_ids ? cell_ids : null ]
        }

    if (!toBool(params.split_by_celltype)) {
        MERGE_BAMS(ch_for_merge.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta, bams, bais ] })

        ch_raw_bams = MERGE_BAMS.out.bam
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_for_merge.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta.id, annot, called_bc, cell_ids ] },
                by: 0
            )
            .map { id, meta, bam, bai, annot, called_bc, cell_ids -> [ meta, bam, bai, annot, called_bc, cell_ids ] }
    } else {
        // ── Cell type as sample: split each lane, then merge ─────────────
        // 10x barcodes are reused across samples, so each lane BAM is split by
        // its own annotation before merging per (sample_id, celltype); a merged
        // BAM is never split by barcode. meta.id = lane_id keeps sinto outputs
        // of lanes of one sample apart.
        ch_lane_for_ct = ch_lane_bams.map { meta, bam, bai, annot, called_bc, cell_id ->
            [ meta + [ id: (meta.lane_id ?: meta.id), arm: meta.id ], bam, bai, annot ]
        }
        SINTO_CELLTYPE_PERLANE(ch_lane_for_ct)

        ch_ct_frags = SINTO_CELLTYPE_PERLANE.out.celltype_bams.flatMap { meta, bam_list, bai_list ->
            def bams = bam_list instanceof List ? bam_list : [bam_list]
            def bais = bai_list instanceof List ? bai_list : [bai_list]
            [bams, bais].transpose().collect { b, bi ->
                def bn = b.baseName.replaceAll('\\.celltype$', '')
                def ct = bn.replaceFirst('^' + java.util.regex.Pattern.quote(meta.id) + '_', '')
                [ [meta.arm, ct], b, bi, (meta.platform ?: 'unknown') ]
            }
        }

        ch_ct_grouped = ch_ct_frags
            .groupTuple(by: 0)
            .map { key, bams, bais, plats ->
                def arm = key[0]; def ct = key[1]
                def ct_meta = [ id: "${arm}__${ct}", arm: arm, celltype: ct, platform: plats[0] ]
                [ ct_meta, bams, bais ]
            }
        MERGE_CELLTYPE_BAMS(ch_ct_grouped)

        ch_raw_bams = MERGE_CELLTYPE_BAMS.out.bam.map { meta, bam, bai ->
            [ meta, bam, bai, file('NO_FILE'), file('NO_FILE'), null ]
        }
    }

    // Per-cell BAM platforms (Smart-seq2) use SELECTIVE_MERGE; barcode platforms use sinto.
    ch_percell_samples = ch_for_merge.filter { meta, bams, bais, annot, called_bc, cell_ids -> cell_ids != null }
    ch_barcode_samples = ch_for_merge.filter { meta, bams, bais, annot, called_bc, cell_ids -> cell_ids == null }

    // ── Step 2: Infer strandedness ──────────────────────────────────────
    refseq_bed = params.refseq_bed ? file(params.refseq_bed) : null
    INFER_STRANDEDNESS(
        ch_raw_bams.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta, bam, bai ] },
        refseq_bed
    )
    ch_strandedness = INFER_STRANDEDNESS.out.strandedness

    // ── Step 3: Generate all BAM type variants ──────────────────────────
    // Labeled channels [meta, bam, bai, bam_label]:
    //                    none        coordinate (GATK)   umi (umi_tools)
    //   sample:          sample      sampleRmDup         sampleUmiDup
    //   goodCells:       goodCells   goodCellsRmDup      goodCellsUmiDup
    //   annotated:       annotated   annotatedRmDup      annotatedUmiDup
    //   celltype:        celltype    celltypeRmDup       celltypeUmiDup
    //   recal / recalq20: coordinate dedup + BQSR, called at params.qual / 20
    //   recaldefault:    recal BAM with GATK defaults (GATK only)
    // A tier is built only if it is in --bam_types or is --percell_bam_type.

    def requested_bam_types = params.bam_types.tokenize(',')
    def requested_callers   = params.callers.tokenize(',')
    def requested_dedup_strategies = params.dedup.tokenize(',')

    def known_bam_types = [
        'sample', 'sampleRmDup', 'sampleUmiDup', 'recal', 'recalq20',
        'goodCells', 'goodCellsRmDup', 'goodCellsUmiDup',
        'annotated', 'annotatedRmDup', 'annotatedUmiDup',
        'celltype', 'celltypeRmDup', 'celltypeUmiDup'
    ]
    def unknown_bam_types = requested_bam_types.findAll { !(it in known_bam_types) }
    if (unknown_bam_types) {
        error "--bam_types has unknown value(s): ${unknown_bam_types.join(', ')}. Allowed: ${known_bam_types.join(', ')}"
    }
    def known_callers = ['bcftools', 'reditools3', 'gatk', 'jacusa2', 'scallele', 'cellsnp']
    def unknown_callers = requested_callers.findAll { !(it in known_callers) }
    if (unknown_callers) {
        error "--callers has unknown value(s): ${unknown_callers.join(', ')}. Allowed: ${known_callers.join(', ')}" +
              (('giremi' in unknown_callers) ? ". GIREMI was removed (human/pm/METHOD_DECISIONS.md)." : '')
    }
    def pc_bam_type_req = toBool(params.percell_calling) ? params.percell_bam_type : null
    def wants = { String t -> (t in requested_bam_types) || t == pc_bam_type_req }
    def wanted_types = known_bam_types.findAll { wants.call(it) }
    def need_coord = wanted_types.findAll { it.endsWith('RmDup') || it.startsWith('recal') }
    if (need_coord && !('coordinate' in requested_dedup_strategies)) {
        error "bam_type(s) ${need_coord.join(', ')} need --dedup to include 'coordinate'"
    }
    def need_umi = wanted_types.findAll { it.endsWith('UmiDup') }
    if (need_umi && !('umi' in requested_dedup_strategies)) {
        error "bam_type(s) ${need_umi.join(', ')} need --dedup to include 'umi'"
    }
    def wants_goodcells = wanted_types.any { it.startsWith('goodCells') }
    def tier_types = wanted_types.findAll { it.startsWith('goodCells') || it.startsWith('annotated') || it.startsWith('celltype') }
    if (toBool(params.split_by_celltype) && tier_types) {
        log.warn "--split_by_celltype: each input is already one cell type, so ${tier_types.join(', ')} are not built. Use sample-level bam_types."
    }

    ch_base_bams = ch_raw_bams.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta, bam, bai ] }

    // goodCells = annotated ∪ QC-passing unannotated
    if (toBool(params.remap)) {
        ch_remap_mats = STARSOLO_ALIGN.out.solo_out
            .map { meta, solo -> [ meta.id, findSoloFiltered(solo) ] }
            .groupTuple()
            .map { id, mats -> [ id, mats.findAll { it }.collect { it.toString() }.unique().join(',') ] }
        ch_qc_input = ch_qc_from_sheet
            .map { meta, annot, called, qc, matrix_str, cells ->
                [ meta.id, meta, annot, called, qc, cells ]
            }
            .combine(ch_remap_mats, by: 0)
            .map { id, meta, annot, called, qc, cells, matrix_str ->
                [ meta, annot, called, qc, matrix_str, cells ]
            }
    } else {
        ch_qc_input = ch_qc_from_sheet
    }
    def run_goodcells_qc = wants_goodcells && !toBool(params.split_by_celltype)
    COMPUTE_GOOD_CELLS(run_goodcells_qc ? ch_qc_input : Channel.empty(), genome_gtf)
    ch_good_bc = COMPUTE_GOOD_CELLS.out.barcodes
        .map { meta, tsv -> [ meta.id, tsv ] }

    // ── 3a. sample — raw BAMs (no dedup) ────────────────────────────────
    ch_sample = ch_base_bams.map { meta, bam, bai -> [ meta, bam, bai, 'sample' ] }

    ch_percell_src = toBool(params.split_by_celltype) ? Channel.empty() : ch_percell_samples

    // ── 3b. Coordinate-deduped BAMs (GATK MarkDuplicates) ──────────────
    //   Per-cell platforms dedup each cell, then merge: SS2 reads carry no cell
    //   tag, so a merged-BAM dedup would collapse fragments across cells.
    // 3b and 3c run per chromosome and per GEM well (@RG SM). Samples pooling GEM
    // wells without @RG take a prebuilt BAM (rmdup_bam / umidup_bam), read as a
    // side channel so `meta`, and every cached task hash, is unchanged.
    def dedup_chroms = params.main_chroms
        ? file(params.main_chroms).readLines().collect { it.trim() }.findAll { it } + ['_other']
        : ['_all']
    def dedup_main = params.main_chroms ? file(params.main_chroms) : file('NO_FILE')
    def route_prebuilt = { ch_bams, column ->
        def ch_pre = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row[column] }
            .map { row -> [ row.sample_id, file(row[column], checkIfExists: true) ] }
            .unique { it[0] }
        ch_bams
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .join(ch_pre, remainder: true)
            .filter { id, meta, bam, bai, pre -> meta != null }
            .branch { id, meta, bam, bai, pre ->
                prebuilt: pre != null
                    return [ meta, pre, file("${pre}.bai", checkIfExists: true) ]
                scatter: true
                    return [ meta, bam, bai ]
            }
    }
    // Emits a sample as soon as all of its dedup_chroms pieces are done.
    def scatter_gather_key = { ch_out -> ch_out
        .map { row -> [ groupKey(row[0], dedup_chroms.size()) ] + row[2..-1] }
        .groupTuple()
        .map { row -> [ row[0].getGroupTarget() ] + row[1..-1] }
    }

    if ('coordinate' in requested_dedup_strategies) {
        ch_coord_route = route_prebuilt.call(ch_raw_bams
            .filter { meta, bam, bai, annot, called_bc, cell_ids -> cell_ids == null }
            .map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta, bam, bai ] }, 'rmdup_bam')
        MARKDUP_CHROM(ch_coord_route.scatter.combine(Channel.fromList(dedup_chroms)), dedup_main)
        MARKDUP_GATHER(scatter_gather_key.call(MARKDUP_CHROM.out.bam))
        MARKDUP_PERCELL_SS2(ch_percell_src.map { meta, bams, bais, annot, called_bc, cell_ids ->
            [ meta, bams, bais, cell_ids ]
        })
        ch_coord_dedup = MARKDUP_GATHER.out.bam.mix(ch_coord_route.prebuilt).mix(MARKDUP_PERCELL_SS2.out.merged)
        ch_percell_rmdup = MARKDUP_PERCELL_SS2.out.bams
            .map { meta, bams, bais, cell_ids ->
                def b  = (bams instanceof List ? bams : [bams]).sort(false) { f -> f.name }
                def bi = (bais instanceof List ? bais : [bais]).sort(false) { f -> f.name }
                [ meta.id, meta, b, bi, cell_ids ]
            }
            .combine(ch_percell_src.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta.id, annot ] }, by: 0)
    }

    // ── 3c. UMI-deduped BAMs (umi_tools) ───────────────────────────────
    if ('umi' in requested_dedup_strategies) {
        ch_umi_route = route_prebuilt.call(ch_base_bams, 'umidup_bam')
        UMI_DEDUP_CHROM(ch_umi_route.scatter.combine(Channel.fromList(dedup_chroms)), dedup_main)
        UMI_DEDUP_GATHER(scatter_gather_key.call(UMI_DEDUP_CHROM.out.bam))
        ch_umi_dedup = UMI_DEDUP_GATHER.out.bam.mix(ch_umi_route.prebuilt)
    }

    ch_sampleRmDup  = ('coordinate' in requested_dedup_strategies)
                ? ch_coord_dedup.map { meta, bam, bai -> [ meta, bam, bai, 'sampleRmDup' ] }
                : Channel.empty()
    ch_sampleUmiDup = ('umi' in requested_dedup_strategies)
                ? ch_umi_dedup.map { meta, bam, bai -> [ meta, bam, bai, 'sampleUmiDup' ] }
                : Channel.empty()

    // ── 3d. Recalibrated BAMs (coordinate dedup + SplitNCigar + BQSR) ──
    // main_chroms keeps GATK off extra / renamed scaffolds, which abort it.
    main_chroms = params.main_chroms ? file(params.main_chroms) : file('NO_FILE')
    def needs_recal = ('recal' in requested_bam_types || 'recalq20' in requested_bam_types)
    if (needs_recal && 'coordinate' in requested_dedup_strategies) {
        GATK_SPLITNCIGAR(ch_coord_dedup, genome_fasta, main_chroms)
        GATK_BQSR(
            GATK_SPLITNCIGAR.out.bam,
            genome_fasta,
            params.known_sites_vcf ? file(params.known_sites_vcf) : [],
            params.known_sites_vcf ? file(params.known_sites_vcf + '.tbi') : [],
            main_chroms
        )
        ch_recal = GATK_BQSR.out.bam.map { meta, bam, bai -> [ meta, bam, bai, 'recal' ] }
        ch_recalq20 = GATK_BQSR.out.bam.map { meta, bam, bai ->
            [ meta + [qual_override: 20], bam, bai, 'recalq20' ]
        }
        ch_gatk_default = GATK_BQSR.out.bam.map { meta, bam, bai -> [ meta, bam, bai, 'recaldefault' ] }
    } else {
        ch_recal = Channel.empty()
        ch_recalq20 = Channel.empty()
        ch_gatk_default = Channel.empty()
    }

    // ── 3e. goodCells ───────────────────────────────────────────────────
    // With --split_by_celltype each input is already one cell type: emptying
    // ch_raw_barcode and ch_percell_src disables 3e-3g.
    ch_raw_barcode = toBool(params.split_by_celltype)
        ? Channel.empty()
        : ch_raw_bams.filter { meta, bam, bai, annot, called_bc, cell_ids -> cell_ids == null }

    // An empty input channel wires a process without launching a task.
    def gate = { List labels, ch -> labels.any { wants.call(it) } ? ch : Channel.empty() }

    // Build the cell-type tier only when asked: celltype_split.nf errors on
    // barcodes assigned to >1 cell type, which merged multi-sample BAMs always
    // have (10x barcodes are reused across samples).
    def wants_celltype = ['celltype', 'celltypeRmDup', 'celltypeUmiDup'].any { wants.call(it) }
    ch_ct_source = wants_celltype ? ch_raw_barcode : Channel.empty()
    if (!wants_celltype) {
        log.info "No celltype bam_type requested — skipping the cell-type split tier."
    }

    ch_for_called_sinto = ch_raw_barcode
        .map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id, meta, bam, bai ] }
        .combine(ch_good_bc, by: 0)
        .map { id, meta, bam, bai, bc -> [ meta, bam, bai, bc, 'goodCells' ] }
    SINTO_GOOD_CELLS(gate.call(['goodCells'], ch_for_called_sinto))
    ch_goodCells_sinto = SINTO_GOOD_CELLS.out.bam.map { meta, bam, bai ->
        [ meta, bam, bai, 'goodCells' ]
    }

    ch_for_called_selmerge = ch_percell_src
        .map { meta, bams, bais, annot, called_bc, cell_ids ->
            [ meta.id, meta, bams, bais, cell_ids ]
        }
        .combine(ch_good_bc, by: 0)
        .map { id, meta, bams, bais, cell_ids, bc ->
            [ meta, bams, bais, cell_ids, bc, 'goodCells', 'filter' ]
        }
    SELECTIVE_MERGE(gate.call(['goodCells'], ch_for_called_selmerge))
    ch_goodCells_selmerge = SELECTIVE_MERGE.out.bam.map { meta, bam, bai ->
        [ meta, bam instanceof List ? bam[0] : bam, bai instanceof List ? bai[0] : bai, 'goodCells' ]
    }

    ch_goodCells = ch_goodCells_sinto.mix(ch_goodCells_selmerge)

    // goodCells + coordinate dedup
    if ('coordinate' in requested_dedup_strategies) {
        ch_coord_with_called = ch_coord_dedup
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_raw_barcode.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id ] },
                by: 0
            )
            .combine(ch_good_bc, by: 0)
            .map { id, meta, bam, bai, bc -> [ meta, bam, bai, bc, 'goodCellsRmDup' ] }
        SINTO_GOOD_CELLS_RMDUP(gate.call(['goodCellsRmDup'], ch_coord_with_called))
        SELECTIVE_MERGE_RMDUP(gate.call(['goodCellsRmDup'],
            ch_percell_rmdup
                .map { id, meta, bams, bais, cell_ids, annot -> [ id, meta, bams, bais, cell_ids ] }
                .combine(ch_good_bc, by: 0)
                .map { id, meta, bams, bais, cell_ids, bc ->
                    [ meta, bams, bais, cell_ids, bc, 'goodCellsRmDup', 'filter' ]
                }))
        ch_goodCellsRmDup = SINTO_GOOD_CELLS_RMDUP.out.bam
            .map { meta, bam, bai -> [ meta, bam, bai, 'goodCellsRmDup' ] }
            .mix(SELECTIVE_MERGE_RMDUP.out.bam.map { meta, bam, bai ->
                [ meta, bam instanceof List ? bam[0] : bam, bai instanceof List ? bai[0] : bai, 'goodCellsRmDup' ]
            })
    } else {
        ch_goodCellsRmDup = Channel.empty()
    }

    // goodCells + UMI dedup (barcode platforms only; per-cell platforms have no UMIs)
    if ('umi' in requested_dedup_strategies) {
        ch_umi_with_called = ch_umi_dedup
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_raw_barcode.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id ] },
                by: 0
            )
            .combine(ch_good_bc, by: 0)
            .map { id, meta, bam, bai, bc -> [ meta, bam, bai, bc, 'goodCellsUmiDup' ] }
        SINTO_GOOD_CELLS_UMIDUP(gate.call(['goodCellsUmiDup'], ch_umi_with_called))
        ch_goodCellsUmiDup = SINTO_GOOD_CELLS_UMIDUP.out.bam.map { meta, bam, bai ->
            [ meta, bam, bai, 'goodCellsUmiDup' ]
        }
    } else {
        ch_goodCellsUmiDup = Channel.empty()
    }

    // ── 3f. annotated — filtered to annotation barcodes ────────────────
    ch_for_annot_sinto = ch_raw_barcode.map { meta, bam, bai, annot, called_bc, cell_ids ->
        [ meta, bam, bai, annot, 'annotated' ]
    }
    SINTO_ANNOTATED(gate.call(['annotated'], ch_for_annot_sinto))
    ch_annotated_sinto = SINTO_ANNOTATED.out.bam.map { meta, bam, bai ->
        [ meta, bam, bai, 'annotated' ]
    }

    ch_for_annot_selmerge = ch_percell_src.map { meta, bams, bais, annot, called_bc, cell_ids ->
        [ meta, bams, bais, cell_ids, annot, 'annotated', 'filter' ]
    }
    SELECTIVE_MERGE_ANNOT(gate.call(['annotated'], ch_for_annot_selmerge))
    ch_annotated_selmerge = SELECTIVE_MERGE_ANNOT.out.bam.map { meta, bam, bai ->
        [ meta, bam instanceof List ? bam[0] : bam, bai instanceof List ? bai[0] : bai, 'annotated' ]
    }

    ch_annotated = ch_annotated_sinto.mix(ch_annotated_selmerge)

    // annotated + coordinate dedup
    if ('coordinate' in requested_dedup_strategies) {
        ch_coord_with_annot = ch_coord_dedup
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_raw_barcode.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id, annot ] },
                by: 0
            )
            .map { id, meta, bam, bai, annot -> [ meta, bam, bai, annot, 'annotatedRmDup' ] }
        SINTO_ANNOTATED_RMDUP(gate.call(['annotatedRmDup'], ch_coord_with_annot))
        SELECTIVE_MERGE_ANNOT_RMDUP(gate.call(['annotatedRmDup'],
            ch_percell_rmdup.map { id, meta, bams, bais, cell_ids, annot ->
                [ meta, bams, bais, cell_ids, annot, 'annotatedRmDup', 'filter' ]
            }))
        ch_annotatedRmDup = SINTO_ANNOTATED_RMDUP.out.bam
            .map { meta, bam, bai -> [ meta, bam, bai, 'annotatedRmDup' ] }
            .mix(SELECTIVE_MERGE_ANNOT_RMDUP.out.bam.map { meta, bam, bai ->
                [ meta, bam instanceof List ? bam[0] : bam, bai instanceof List ? bai[0] : bai, 'annotatedRmDup' ]
            })
    } else {
        ch_annotatedRmDup = Channel.empty()
    }

    // annotated + UMI dedup
    if ('umi' in requested_dedup_strategies) {
        ch_umi_with_annot = ch_umi_dedup
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_raw_barcode.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id, annot ] },
                by: 0
            )
            .map { id, meta, bam, bai, annot -> [ meta, bam, bai, annot, 'annotatedUmiDup' ] }
        SINTO_ANNOTATED_UMIDDUP(gate.call(['annotatedUmiDup'], ch_umi_with_annot))
        ch_annotatedUmiDup = SINTO_ANNOTATED_UMIDDUP.out.bam.map { meta, bam, bai ->
            [ meta, bam, bai, 'annotatedUmiDup' ]
        }
    } else {
        ch_annotatedUmiDup = Channel.empty()
    }

    // ── 3g. Cell-type BAMs (one BAM per cell type) ───────────────────────

    SINTO_EXTRACT_CELLTYPE(gate.call(['celltype'], ch_ct_source.map { meta, bam, bai, annot, called_bc, cell_ids ->
        [ meta, bam, bai, annot ]
    }))
    ch_celltype_sinto = flattenCelltype(SINTO_EXTRACT_CELLTYPE.out.celltype_bams, 'celltype')

    ch_for_ct_selmerge = ch_percell_src.map { meta, bams, bais, annot, called_bc, cell_ids ->
        [ meta, bams, bais, cell_ids, annot, 'celltype', 'celltype' ]
    }
    SELECTIVE_MERGE_CT(gate.call(['celltype'], ch_for_ct_selmerge))
    ch_celltype_selmerge = flattenSelectiveCelltype(SELECTIVE_MERGE_CT.out.bam, 'celltype')

    ch_celltype = ch_celltype_sinto.mix(ch_celltype_selmerge)

    // celltype + coordinate dedup
    if ('coordinate' in requested_dedup_strategies) {
        ch_coord_for_ct = ch_coord_dedup
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_ct_source.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id, annot ] },
                by: 0
            )
            .map { id, meta, bam, bai, annot -> [ meta, bam, bai, annot ] }
        SINTO_CELLTYPE_RMDUP(gate.call(['celltypeRmDup'], ch_coord_for_ct))
        SELECTIVE_MERGE_CT_RMDUP(gate.call(['celltypeRmDup'],
            ch_percell_rmdup.map { id, meta, bams, bais, cell_ids, annot ->
                [ meta, bams, bais, cell_ids, annot, 'celltypeRmDup', 'celltype' ]
            }))
        ch_celltypeRmDup = flattenCelltype(SINTO_CELLTYPE_RMDUP.out.celltype_bams, 'celltypeRmDup')
            .mix(flattenSelectiveCelltype(SELECTIVE_MERGE_CT_RMDUP.out.bam, 'celltypeRmDup'))
    } else {
        ch_celltypeRmDup = Channel.empty()
    }

    // celltype + UMI dedup
    if ('umi' in requested_dedup_strategies) {
        ch_umi_for_ct = ch_umi_dedup
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_ct_source.map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id, annot ] },
                by: 0
            )
            .map { id, meta, bam, bai, annot -> [ meta, bam, bai, annot ] }
        SINTO_CELLTYPE_UMIDDUP(gate.call(['celltypeUmiDup'], ch_umi_for_ct))
        ch_celltypeUmiDup = flattenCelltype(SINTO_CELLTYPE_UMIDDUP.out.celltype_bams, 'celltypeUmiDup')
    } else {
        ch_celltypeUmiDup = Channel.empty()
    }

    // ── Assemble BAM channels for variant calling ───────────────────────
    ch_all_labeled_bams = Channel.empty()
    if ('sample'             in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_sample)
    if ('sampleRmDup'              in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_sampleRmDup)
    if ('sampleUmiDup'             in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_sampleUmiDup)
    if ('recal'              in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_recal)
    if ('recalq20'           in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_recalq20)
    if ('goodCells'        in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_goodCells)
    if ('goodCellsRmDup'   in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_goodCellsRmDup)
    if ('goodCellsUmiDup'  in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_goodCellsUmiDup)
    if ('annotated'          in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_annotated)
    if ('annotatedRmDup'     in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_annotatedRmDup)
    if ('annotatedUmiDup'    in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_annotatedUmiDup)
    if ('celltype'           in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_celltype)
    if ('celltypeRmDup'      in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_celltypeRmDup)
    if ('celltypeUmiDup'     in requested_bam_types) ch_all_labeled_bams = ch_all_labeled_bams.mix(ch_celltypeUmiDup)

    // ── Step 4: Variant calling ─────────────────────────────────────────
    ch_bams_with_strand = ch_all_labeled_bams
        .map { meta, bam, bai, label -> [ meta.parent_id ?: meta.id, meta, bam, bai, label ] }
        .combine(
            ch_strandedness.map { meta, strand -> [ meta.id, strand ] },
            by: 0
        )
        .map { id, meta, bam, bai, label, strand -> [ meta, bam, bai, label, strand ] }

    main_chroms   = params.main_chroms   ? file(params.main_chroms) : file('NO_FILE')

    if ('bcftools' in requested_callers) {
        BCFTOOLS_MPILEUP(ch_all_labeled_bams, genome_fasta, main_chroms)
    }
    if ('reditools3' in requested_callers) {
        REDITOOLS3(ch_bams_with_strand, genome_fasta, main_chroms)
    }
    if ('gatk' in requested_callers) {
        // recaldefault (GATK defaults on the recal BAM) is run by GATK only.
        GATK_HAPLOTYPECALLER(ch_all_labeled_bams.mix(ch_gatk_default), genome_fasta, main_chroms)
    }
    def sc_chroms = params.main_chroms
        ? file(params.main_chroms).readLines().collect { it.trim() }.findAll { it }
        : file("${params.genome_fasta}.fai", checkIfExists: true).readLines().collect { it.tokenize('\t')[0] }
    def ch_jacusa2_variants = Channel.empty()
    if ('jacusa2' in requested_callers) {
        if (toBool(params.jacusa2_scatter)) {
            if (!params.main_chroms)
                log.warn "JACUSA2: no --main_chroms, scattering over all ${sc_chroms.size()} reference contigs"
            JACUSA2_CHROM(ch_bams_with_strand.combine(Channel.fromList(sc_chroms)), genome_fasta)
            JACUSA2_GATHER(
                JACUSA2_CHROM.out.chrom_out
                    .map { meta, label, chrom, out -> [ groupKey([meta, label], sc_chroms.size()), chrom, out ] }
                    .groupTuple()
                    .map { key, chroms, outs -> def k = key.getGroupTarget(); [ k[0], k[1], chroms, outs ] },
                sc_chroms)
            ch_jacusa2_variants = JACUSA2_GATHER.out.variants
        } else {
            JACUSA2_CALL1(ch_bams_with_strand, genome_fasta, main_chroms)
            ch_jacusa2_variants = JACUSA2_CALL1.out.variants
        }
    }
    if ('cellsnp' in requested_callers) {
        CELLSNP_LITE(ch_all_labeled_bams, genome_fasta, main_chroms)
    }
    if ('scallele' in requested_callers) {
        if (!params.main_chroms)
            log.warn "scAllele: no --main_chroms, scattering over all ${sc_chroms.size()} reference contigs"
        SCALLELE_CHROM(ch_bams_with_strand.combine(Channel.fromList(sc_chroms)), genome_fasta)
        SCALLELE_GATHER(
            SCALLELE_CHROM.out.chrom_out
                .map { meta, label, chrom, snv, vcf, rc -> [ groupKey([meta, label], sc_chroms.size()), chrom, snv, vcf, rc ] }
                .groupTuple()
                .map { key, chroms, snvs, vcfs, rcs -> def k = key.getGroupTarget(); [ k[0], k[1], chroms, snvs, vcfs, rcs ] },
            sc_chroms)
    }

    // ── Step 5: Downstream R processing ─────────────────────────────────
    snp_file      = params.snp_file     ? file(params.snp_file) : file('NO_FILE')
    strand_annots = params.strand_annots ? file(params.strand_annots) : file('NO_FILE')
    alu_annots    = params.alu_annots    ? file(params.alu_annots) : file('NO_FILE')

    // 5th element: per-sample genetic-variant blacklist (bcftools / gatk 1/1 calls;
    // NO_FILE otherwise). scAllele's PASS is not a genotype call: well-supported
    // edits PASS too.
    no_blacklist = file('NO_FILE')
    def join_blacklist = { var_ch, blk_ch ->
        var_ch.map { m, c, bl, v -> [ [m.id, bl], m, c, bl, v ] }
              .join( blk_ch.map { m, c, bl, blk -> [ [m.id, bl], blk ] } )
              .map { key, m, c, bl, v, blk -> [ m, c, bl, v, blk ] }
    }
    ch_all_variants = Channel.empty()
    if ('bcftools'    in requested_callers) ch_all_variants = ch_all_variants.mix(join_blacklist.call(BCFTOOLS_MPILEUP.out.variants, BCFTOOLS_MPILEUP.out.snp_blacklist))
    if ('reditools3'  in requested_callers) ch_all_variants = ch_all_variants.mix(REDITOOLS3.out.variants.map { m, c, bl, v -> [ m, c, bl, v, no_blacklist ] })
    if ('gatk'        in requested_callers) ch_all_variants = ch_all_variants.mix(join_blacklist.call(GATK_HAPLOTYPECALLER.out.variants, GATK_HAPLOTYPECALLER.out.snp_blacklist))
    if ('jacusa2'     in requested_callers) ch_all_variants = ch_all_variants.mix(ch_jacusa2_variants.map { m, c, bl, v -> [ m, c, bl, v, no_blacklist ] })
    if ('scallele'    in requested_callers) ch_all_variants = ch_all_variants.mix(SCALLELE_GATHER.out.variants.map { m, c, bl, v -> [ m, c, bl, v, no_blacklist ] })
    if ('cellsnp'     in requested_callers) ch_all_variants = ch_all_variants.mix(CELLSNP_LITE.out.variants.map { m, c, bl, v -> [ m, c, bl, v, no_blacklist ] })
    // ── Benchmark axis: per-sample SNP subtraction ('on' needs a blacklist) ──
    def snp_modes = params.snp_subtract_modes.tokenize(',')
    ch_variants_snpmodes = ch_all_variants.flatMap { meta, caller, bam_label, variants, blacklist ->
        def rows = []
        if ('off' in snp_modes) rows << [ meta, caller, bam_label, variants, 'off', no_blacklist ]
        if ('on'  in snp_modes && blacklist.name != 'NO_FILE') rows << [ meta, caller, bam_label, variants, 'on', blacklist ]
        return rows
    }

    def requested_strand_modes = params.strand_modes.tokenize(',')
    ch_strand_modes = Channel.fromList(requested_strand_modes)

    // 'hybrid' needs I16/SB read-strand counts (bcftools, gatk); for other callers it
    // would duplicate 'annotation'. 'caller' needs a caller-emitted strand
    // (reditools3, jacusa2, scallele); bcftools and gatk report "+" only.
    // 'caller' is deliberately not gated on library strandedness: on unstranded
    // libraries it shows what users of that tool would get.
    def hybrid_callers        = ['bcftools', 'gatk']
    def caller_strand_callers = ['reditools3', 'jacusa2', 'scallele']
    ch_variants_x_strand = ch_variants_snpmodes
        .combine(ch_strand_modes)
        .filter { meta, caller, bam_label, variants, snp_sub, blacklist, strand_mode ->
            strand_mode == 'annotation' ||
            (strand_mode == 'hybrid' && caller in hybrid_callers) ||
            (strand_mode == 'caller' && caller in caller_strand_callers)
        }
        .map { meta, caller, bam_label, variants, snp_sub, blacklist, strand_mode ->
            [ meta, caller, bam_label, variants, strand_mode, blacklist, snp_sub ]
        }

    // scAllele has its own copy of DOWNSTREAM_R.
    ch_ds_split = ch_variants_x_strand.branch { row ->
        scallele: row[1] == 'scallele'
        other:    true
    }
    DOWNSTREAM_R(
        ch_ds_split.other,
        genome_fasta,
        genome_gtf,
        snp_file,
        strand_annots,
        alu_annots,
        main_chroms
    )
    DOWNSTREAM_R_SCALLELE(
        ch_ds_split.scallele,
        genome_fasta,
        genome_gtf,
        snp_file,
        strand_annots,
        alu_annots,
        main_chroms
    )
    ch_ds_filtered = DOWNSTREAM_R.out.filtered.mix(DOWNSTREAM_R_SCALLELE.out.filtered)
    ch_ds_pre100   = DOWNSTREAM_R.out.pre100.mix(DOWNSTREAM_R_SCALLELE.out.pre100)
    ch_ds_stepwise = DOWNSTREAM_R.out.stepwise.mix(DOWNSTREAM_R_SCALLELE.out.stepwise)
    ch_ds_metrics  = DOWNSTREAM_R.out.metrics.mix(DOWNSTREAM_R_SCALLELE.out.metrics)

    // ── Step 5b: Re-quantify editing level + coverage (strand-aware pileup) ──
    // Uses the bam_type BAM that produced each call. REQUANT_PILEUP, not
    // DOWNSTREAM_R, then publishes filtered / pre100 to 04_filtered.
    if (toBool(params.requant_editing)) {
        ch_bam_lut = ch_all_labeled_bams.mix(ch_gatk_default)
            .map { m, bam, bai, label -> [ [m.id, label], bam, bai ] }

        ch_fp = ch_ds_filtered
            .map { m, c, bl, sm, ss, f -> [ [m.id, c, bl, sm, ss], m, c, bl, sm, ss, f ] }
            .join(
                ch_ds_pre100.map { m, c, bl, sm, ss, p -> [ [m.id, c, bl, sm, ss], p ] }
            )
            .map { key, m, c, bl, sm, ss, f, p -> [ [m.id, bl], m, c, bl, sm, ss, f, p ] }

        ch_requant_in = ch_fp
            .combine(ch_bam_lut, by: 0)
            .map { key, m, c, bl, sm, ss, f, p, bam, bai -> [ m.parent_id ?: m.id, m, c, bl, sm, ss, f, p, bam, bai ] }
            .combine(ch_strandedness.map { m, s -> [ m.id, s ] }, by: 0)
            .map { id, m, c, bl, sm, ss, f, p, bam, bai, s -> [ m, c, bl, sm, ss, f, p, bam, bai, s ] }

        REQUANT_PILEUP(
            ch_requant_in,
            genome_fasta,
            file("${params.genome_fasta}.fai")
        )
    }

    ch_ds_stepwise
        .map { meta, caller, bam_label, strand_mode, tsv -> tsv }
        .collectFile(
            name: 'all_stepwise_counts.tsv',
            storeDir: "${params.outdir}/summary",
            keepHeader: true,
            skip: 1
        )

    ch_all_metrics = ch_ds_metrics
        .map { meta, caller, bam_label, strand_mode, tsv -> tsv }
        .collectFile(
            name: 'all_run_metrics.tsv',
            storeDir: "${params.outdir}/summary",
            keepHeader: true,
            skip: 1
        )

    COLLATE_RESULTS(ch_all_metrics)

    // ── REDIportal overlap summary (human) ──────────────────────────────────
    if (params.rediportal) {
        ch_pre100_all = ch_ds_pre100
            .map { m, c, bl, sm, ss, p -> p }
            .collect()
        REDIPORTAL_OVERLAP(ch_pre100_all, file(params.rediportal))
    }

    ch_coverage = Channel.empty()
    ch_coverage = ch_coverage.mix(SINTO_GOOD_CELLS.out.coverage)
    ch_coverage = ch_coverage.mix(SINTO_ANNOTATED.out.coverage)
    if ('coordinate' in requested_dedup_strategies) {
        ch_coverage = ch_coverage.mix(SINTO_GOOD_CELLS_RMDUP.out.coverage)
        ch_coverage = ch_coverage.mix(SINTO_ANNOTATED_RMDUP.out.coverage)
    }
    if ('umi' in requested_dedup_strategies) {
        ch_coverage = ch_coverage.mix(SINTO_GOOD_CELLS_UMIDUP.out.coverage)
        ch_coverage = ch_coverage.mix(SINTO_ANNOTATED_UMIDDUP.out.coverage)
    }
    ch_coverage
        .map { meta, tsv -> tsv }
        .collectFile(
            name: 'barcode_coverage_summary.tsv',
            storeDir: "${params.outdir}/summary",
            keepHeader: true,
            skip: 1
        )

    COMPUTE_GOOD_CELLS.out.summary
        .map { meta, tsv -> tsv }
        .collectFile(
            name: 'barcode_qc_summary.tsv',
            storeDir: "${params.outdir}/summary",
            keepHeader: true,
            skip: 1
        )

    // ── Step 6: Per-cell analyses (barcode platforms only: needs CB tags) ──
    if (toBool(params.percell_analyses)) {
    ch_barcode_base_bams = ch_base_bams
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .combine(
            ch_barcode_samples.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta.id, called_bc ?: file('NO_FILE') ] },
            by: 0
        )
        .map { id, meta, bam, bai, called_bc -> [ meta, bam, bai, called_bc ] }

    // Without called barcodes sinto opens one file per background barcode (> 250k in 10x BAMs).
    SINTO_SPLIT_PERCELL(ch_barcode_base_bams)
    PERCELL_MPILEUP(
        SINTO_SPLIT_PERCELL.out.percell_bams,
        genome_fasta
    )

    if (params.repeat_annots) {
        PERCELL_AEI(
            SINTO_SPLIT_PERCELL.out.percell_bams,
            file(params.repeat_annots),
            genome_fasta
        )
    }
    }  // end params.percell_analyses

    // ── Step 7: Per-cell discovery (optional, second-pass mode) ─────────
    if (toBool(params.percell_calling)) {
        def pc_bam_type = params.percell_bam_type

        // Cell-type tiers are excluded: one sample yields several BAMs there,
        // and splitting them per cell just re-derives the annotated tier.
        def pc_sources = [
            sample:          ch_sample,
            sampleRmDup:     ch_sampleRmDup,
            sampleUmiDup:    ch_sampleUmiDup,
            recal:           ch_recal,
            goodCells:       ch_goodCells,
            goodCellsRmDup:  ch_goodCellsRmDup,
            goodCellsUmiDup: ch_goodCellsUmiDup,
            annotated:       ch_annotated,
            annotatedRmDup:  ch_annotatedRmDup,
            annotatedUmiDup: ch_annotatedUmiDup
        ]
        if (!pc_sources.containsKey(pc_bam_type)) {
            error "--percell_bam_type must be one of: ${pc_sources.keySet().join(', ')} (got '${pc_bam_type}')"
        }
        ch_percell_bams = pc_sources[pc_bam_type].map { meta, bam, bai, label -> [ meta, bam, bai ] }

        // Barcode list keeps background barcodes of the sample BAM out.
        if (toBool(params.split_by_celltype)) {
            ch_percell_in = ch_percell_bams.map { meta, bam, bai -> [ meta, bam, bai, file('NO_FILE') ] }
        } else {
            ch_pc_bc = pc_bam_type.startsWith('goodCells')
                ? ch_good_bc
                : ch_for_merge.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta.id, called_bc ?: annot ] }
            ch_percell_in = ch_percell_bams
                .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
                .combine(ch_pc_bc, by: 0)
                .map { id, meta, bam, bai, bc -> [ meta, bam, bai, bc ] }
        }

        // Split once, then call per cell (avoids one whole-BAM scan per cell).
        SINTO_SPLIT_PERCELL_DISC(ch_percell_in)

        PERCELL_CALLING(
            SINTO_SPLIT_PERCELL_DISC.out.percell_bams,
            genome_fasta,
            file("${params.genome_fasta}.fai"),
            file("${params.genome_fasta}".replaceAll('\\.fa(sta)?$', '.dict')),
            snp_file,
            strand_annots,
            file(params.main_chroms ?: 'NO_FILE')
        )
    }
}
