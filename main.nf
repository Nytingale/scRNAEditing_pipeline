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
params.genome_fasta    = null    // Reference genome FASTA
params.genome_gtf      = null    // Gene annotation GTF
params.star_index      = null    // Pre-built STAR index (optional)
params.snp_file        = null    // Known SNP file — RData (GRanges, preferred) or BED (legacy)
params.known_sites_vcf = null    // Known sites VCF/BED for GATK BQSR (null = skip recalibration)
params.strand_annots   = null    // Pre-built strand annotation RData (optional)
params.alu_annots      = null    // Alu annotations RData for human (optional)
params.rediportal      = null    // REDIportal VCF(.gz) of known editing positions (human). When set, emit summary/rediportal_overlap.tsv (n known vs novel + fractions) for figures.
params.repeat_annots   = null    // Repeat annotations for AEI (optional)
params.refseq_bed      = null    // RefSeq BED for strandedness inference
params.main_chroms     = null    // One-per-line file of chromosomes to analyze (null = all from BAM header)
params.remap           = false   // Re-align with STARsolo
params.dedup           = 'none,coordinate'  // Comma-separated dedup strategies to benchmark: 'none', 'coordinate', 'umi'
params.dedup_cell_tag  = 'CR'          // Cell barcode tag for dedup grouping: 'CR' (raw, recommended) or 'CB' (corrected)
params.callers         = 'bcftools,reditools3,gatk,jacusa2'  // add scallele / cellsnp via --callers; not default
params.bam_types       = 'sample,sampleRmDup,sampleUmiDup,recal,recalq20,goodCells,goodCellsRmDup,goodCellsUmiDup,annotated,annotatedRmDup,annotatedUmiDup,celltype,celltypeRmDup,celltypeUmiDup'
// Cell-type-level discovery: when true, each merged sample BAM is split by cell
// type (via celltype_annot) up-front and EACH cell type becomes its own sample
// (id = <sample>__<celltype>). The full sample-level matrix (dedup none/
// coordinate/umi + BQSR recal/recalq20) then runs per cell type -> split-then-
// dedup and BQSR-per-celltype, with unique output paths (no id collision).
// Use sample-level bam_types (sample,sampleRmDup,sampleUmiDup,recal,recalq20).
params.split_by_celltype   = false
params.min_celltype_cells  = 0   // In split_by_celltype mode, drop cell types with fewer cells (0 = keep all)
params.qual            = 30      // Base/mapping quality threshold
params.adaptive_windows = '0,20,40,60,80,100,120,140,160,180,200'
// Pass 1 of the window step (remove sites with another mismatch type within w1 bp). Off since
// 2026-10-01: only the clustered-sites pass (w2, searched over adaptive_windows) is applied.
params.mismatch_window_filter = false
params.min_cov         = 10      // Minimum read coverage at a site
params.strand_modes    = 'annotation'  // Strand assignment modes to benchmark: 'annotation', 'hybrid', or 'annotation,hybrid'
params.strand_ambiguity_threshold = 0.3  // For hybrid mode: read-strand ambiguity cutoff

// Per-cell discovery parameters (second-pass mode)
params.percell_analyses      = true    // Step 6 per-cell split + mpileup/AEI. Set false for caller-comparison-only runs (avoids the memory-heavy SINTO_SPLIT_PERCELL).
params.percell_calling       = false   // Enable per-cell discovery mode (Step 7)
params.percell_aei           = true    // Phase 1 chunked RNAEditingIndexer inside PERCELL_CALLING. Set false for discovery-only (variant-calling) runs.
params.percell_caller        = null    // Required when percell_calling=true: 'bcftools', 'reditools3', 'gatk', 'jacusa2', 'scallele'
params.percell_bam_type      = null    // Required when percell_calling=true: which BAM type to use
params.percell_chunk_size    = 100     // Cells per RNAEditingIndexer chunk
params.rnaeditingindexer_path = null   // Path to RNAEditingIndexer install dir
params.aei_genome            = null    // RNAEditingIndexer --genome (e.g., 'hg38')
params.aei_regions_bed       = null    // Custom repeat regions BED for AEI
params.aei_snps_bed          = null    // Custom SNPs BED for AEI
params.aei_refseq_bed        = null    // Custom RefSeq BED for AEI
params.aei_expression_bed    = null    // Custom gene expression BED for AEI
params.min_alt               = 3       // Minimum alt allele read count (--alt_rule count)
params.min_edit_level        = 0.01    // Minimum editing level
// Alt-read support rule: 'count' (cov_alt >= min_alt) or 'binom' (alt count
// significantly above the sequencing error rate). `count` is the historical
// default and is kept as the default so existing runs are unchanged.
// min_alt=3 is the BINDING sensitivity threshold, not min_edit_level: it implies
// editing >= 3/coverage, i.e. 8.3% at the median depth of 36 reads, while
// min_edit_level=0.01 binds on only 2.58% of sites. See bin/downstream_filter.R.
params.alt_rule              = 'count'
params.seq_error_rate        = 0.001   // per-base error under the null, for 'binom'
params.alt_binom_alpha       = 0.05    // significance threshold, for 'binom'
params.min_alt_floor         = 2       // never accept a single alt read, for 'binom'
params.jacusa2_min_score     = 1.25    // Minimum JACUSA2 quality score
params.jacusa2_scatter       = false   // true: JACUSA2 as one task per chromosome (JACUSA2_CHROM + JACUSA2_GATHER) instead of one whole-genome JACUSA2_CALL1 task
params.scallele_conda        = '/home/yng7/miniforge3/envs/scAllele'
params.scallele_run_mode     = 'Variant_Caller'  // Variant_Caller | Full | Training
params.cellsnp_bin           = '/home/yng7/miniforge3/envs/cellsnp/bin/cellsnp-lite'  // cellSNP-lite 1.2.3 (--callers cellsnp; not default)
params.snp_subtract_modes    = 'off,on'  // Benchmark axis for per-sample SNP filtering. 'off'=dbSNP only; 'on'=dbSNP + caller's own per-sample genetic-variant blacklist. Only genotyping callers (gatk, bcftools) run 'on'; reditools3/jacusa2/scallele are 'off' only (NA).
params.remove_100pct_editing = true    // Remove sites with 100% editing level
params.save_intermediates    = 0       // Save intermediate tables from step N onwards (0=off)
params.requant_editing       = true    // Re-quantify editing level + coverage at each site from the matching bam_type BAM (strand-aware pileup) and ADD canonical columns (edit_level_pileup, cov_pileup, edited_reads, ref_reads) to the published filtered/pre100 tables. Makes editing levels comparable across callers.
params.requant_min_bq        = 30      // requant: min base quality for mpileup
params.requant_min_mq        = 30      // requant: min mapping quality for mpileup
params.requant_max_depth     = 1000000 // requant: samtools mpileup -d
params.mito_chr              = null    // Mitochondrial chromosome override (auto-detect from GTF)
params.qc_min_features       = 200     // goodCells: min genes detected in unannotated cells
params.qc_max_features       = 6000    // goodCells: max genes detected in unannotated cells
params.qc_max_pct_mt         = 20      // goodCells: max mitochondrial percent
params.qc_min_counts         = 500     // goodCells: min UMIs/counts in unannotated cells


// ── Validate inputs (moved to workflow block for Nextflow >=26 compatibility) ──


// Warn about called_barcodes after samplesheet is parsed (checked in workflow)


// ── Include modules ─────────────────────────────────────────────────────────
include { VALIDATE_BAMS        } from './modules/bam_validate'
include { VALIDATE_FASTQ_READS } from './modules/fastq_validate'
include { STARSOLO_ALIGN       } from './modules/starsolo'
include { STAR_GENOMEGENERATE  } from './modules/starsolo'
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


// (Samplesheet parsing moved inside workflow block for Nextflow >=26 DSL2 compatibility)

// Coerce CLI/params booleans: `--flag false` on the command line arrives as the
// STRING "false", which is truthy in Groovy. Normalize to a real Boolean.
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

// Cell-type tier BAMs get their own id (<sample>__<celltype>, same scheme as
// --split_by_celltype) so caller/filter outputs, which are named and joined by
// meta.id, do not overwrite each other across cell types of one sample.
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

// SELECTIVE_MERGE celltype mode emits <celltype>.bam files (one per cell type).
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

    // Validate per-cell discovery parameters
    if (toBool(params.percell_calling)) {
        if (!params.percell_caller)   { error "--percell_calling requires --percell_caller <bcftools|reditools3|gatk|jacusa2|scallele>" }
        if (!params.percell_bam_type) { error "--percell_calling requires --percell_bam_type <sample|sampleRmDup|annotated|...>" }
        // Per-cell AEI (RNAEditingIndexer) is handled by Step 6 (PERCELL_AEI),
        // not by PERCELL_CALLING; no rnaeditingindexer_path requirement here.
    }
    // ── Parse samplesheet ───────────────────────────────────────────────
    // CSV columns:
    //   sample_id, bam, celltype_annot, called_barcodes, platform, fastq_r1, fastq_r2, cell_id
    //
    // Required: sample_id, celltype_annot
    // Optional: called_barcodes — fallback cell list when no STARsolo/CellRanger
    //           count matrix is available (e.g. Smart-seq2 QC-filtered cell_ids).
    // Optional: count_matrix — 10x-style matrix dir (matrix.mtx + features + barcodes).
    //           Auto-discovered from <bam_prefix>_Solo.out/GeneFull/filtered when omitted.
    // Optional: qc_barcodes — precomputed goodCells list; skips metric QC.
    // Optional: cell_id — cell identifier for per-cell BAM platforms (e.g. Smart-seq2).
    //           When provided, the pipeline uses selective samtools merge (instead of sinto)
    //           for goodCells/annotated/celltype BAM types, merging only the per-cell BAMs
    //           whose cell_id matches the annotation file.
    // One of: bam (if starting from BAMs) or fastq_r1+fastq_r2+platform (if remapping)
    //
    // Lane vs sample (IMPORTANT — correctness):
    //   sample_id is the TRUE biological sample and the sample-level ANALYSIS
    //   partition. Multiple rows sharing a sample_id are LANES / technical
    //   replicates of that SAME sample (same cells → SAME barcodes); each row is
    //   processed independently (align/validate/MD/unique-filter) and then MERGED
    //   (samtools merge) before any analysis. Lanes are therefore NEVER treated as
    //   separate samples; lane_id (from the BAM/FASTQ name) is only a pre-merge
    //   join key and never a partition/output key.
    //   DISTINCT biological samples MUST get DISTINCT sample_ids. Do NOT pool
    //   different samples under one sample_id: 10x barcodes are reused across
    //   samples, so a barcode-based split of the naively-merged BAM would MISASSIGN
    //   reads. To deliberately pool cell types ACROSS samples (e.g. teratoma arms
    //   T1-T4), use --split_by_celltype, which splits each row (sample) by cell
    //   type BEFORE merging (collision-safe) rather than splitting a merged BAM.

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

    // Per-sample inputs for goodCells QC (one process per sample_id)
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
    // Each samplesheet row is a "lane" (or replicate). Multiple rows with the same
    // sample_id are processed independently, then merged. Single-row samples are
    // passed through with a symlink (no merge overhead).

    if (toBool(params.remap)) {
        ch_to_remap = ch_input
            .filter { meta, bam, annot, called_bc, r1, r2, cell_id -> r1 != null && r2 != null }
            .map { meta, bam, annot, called_bc, r1, r2, cell_id -> [ meta, annot, called_bc, r1, r2, cell_id ] }

        // Validate FASTQ read identity (barcode vs cDNA) before alignment
        VALIDATE_FASTQ_READS(ch_to_remap.map { meta, annot, called_bc, r1, r2, cell_id -> [ meta, r1, r2 ] })

        // Rejoin validated FASTQs with annotation, pass to STARSOLO_ALIGN
        ch_validated = VALIDATE_FASTQ_READS.out.validated
            .map { meta, bc_fq, cdna_fq -> [ meta.lane_id ?: meta.id, meta, bc_fq, cdna_fq ] }
            .combine(
                ch_to_remap.map { meta, annot, called_bc, r1, r2, cell_id -> [ meta.lane_id ?: meta.id, annot ] },
                by: 0
            )
            .map { lid, meta, bc_fq, cdna_fq, annot -> [ meta, annot, bc_fq, cdna_fq ] }

        STARSOLO_ALIGN(ch_validated, ch_star_index, genome_fasta, genome_gtf)
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

    // Group per-lane BAMs by sample_id and merge
    // ch_lane_bams: [meta, bam, bai, annot, called_bc, cell_id]
    // After grouping: [sample_meta, [bam1, bam2, ...], [bai1, bai2, ...], annot, called_bc, [cell_id1, ...]]
    ch_for_merge = ch_lane_bams
        .map { meta, bam, bai, annot, called_bc, cell_id ->
            [ meta.id, meta, bam, bai, annot, called_bc, cell_id ]
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, bams, bais, annots, called_bcs, cell_ids ->
            def sample_meta = [id: sample_id, platform: metas[0].platform]
            def annot = annots[0]
            // Lanes of one sample must share ONE annotation (same cells). If a
            // sample_id groups rows with DIFFERENT annotation files, only the first
            // is used for barcode-based tiers on the merged BAM (non-celltype path)
            // — a silent drop that usually means distinct samples were pooled under
            // one sample_id. Flag it (use distinct sample_ids, or --split_by_celltype
            // which honours each row's own annotation).
            if (!toBool(params.split_by_celltype)) {
                def distinct_annots = annots.findAll { it != null }.collect { it.toString() }.unique()
                if (distinct_annots.size() > 1) {
                    log.warn "sample_id '${sample_id}' has ${distinct_annots.size()} different celltype_annot files across its rows; only the first is used for barcode-based tiers. Rows sharing a sample_id should be LANES of ONE sample (same annotation). For distinct samples use distinct sample_ids; to pool cell types across samples use --split_by_celltype."
                }
            }
            def called_bc = called_bcs.find { it != null }
            // cell_ids: list of cell identifiers (Smart-seq2) or list of nulls (10x)
            def has_cell_ids = cell_ids.any { it != null }
            [ sample_meta, bams, bais, annot, called_bc, has_cell_ids ? cell_ids : null ]
        }

    if (!toBool(params.split_by_celltype)) {
        // Sample-level merge (all BAMs per sample_id → one merged BAM)
        MERGE_BAMS(ch_for_merge.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta, bams, bais ] })

        ch_raw_bams = MERGE_BAMS.out.bam
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(
                ch_for_merge.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta.id, annot, called_bc, cell_ids ] },
                by: 0
            )
            .map { id, meta, bam, bai, annot, called_bc, cell_ids -> [ meta, bam, bai, annot, called_bc, cell_ids ] }
    } else {
        // ── Cell-type-level discovery — COLLISION-SAFE split-then-merge ──────
        // CRITICAL: 10x cell barcodes are REUSED across samples/lanes, so
        // splitting a MERGED multi-sample BAM by raw CB would pool reads from
        // different cells that happen to share a barcode → misassigned reads.
        // Instead split EACH per-lane BAM (a single sample/lane, where barcodes
        // ARE unique) by cell type using THAT lane's own annotation, then merge the
        // per-cell-type fragments across lanes per (sample_id, celltype). No
        // barcode is ever matched across wells, so misassignment is impossible by
        // construction. Each cell type then becomes its own sample
        // (id = <sample>__<celltype>) and the full sample-level matrix runs per
        // cell type (split-then-dedup + BQSR-per-celltype), with unique output
        // paths. The old 3e/3f/3g secondary splits are neutralized below
        // (ch_raw_barcode emptied); request sample-level bam_types for these runs.
        // Per-lane meta.id is temporarily set to lane_id so sinto output names /
        // publishDir do not collide between lanes of the same sample.
        ch_lane_for_ct = ch_lane_bams.map { meta, bam, bai, annot, called_bc, cell_id ->
            [ meta + [ id: (meta.lane_id ?: meta.id), arm: meta.id ], bam, bai, annot ]
        }
        SINTO_CELLTYPE_PERLANE(ch_lane_for_ct)

        // Flatten per-lane cell-type BAMs; key by (sample_id, celltype) for pooling.
        ch_ct_frags = SINTO_CELLTYPE_PERLANE.out.celltype_bams.flatMap { meta, bam_list, bai_list ->
            def bams = bam_list instanceof List ? bam_list : [bam_list]
            def bais = bai_list instanceof List ? bai_list : [bai_list]
            [bams, bais].transpose().collect { b, bi ->
                def bn = b.baseName.replaceAll('\\.celltype$', '')
                def ct = bn.replaceFirst('^' + java.util.regex.Pattern.quote(meta.id) + '_', '')
                [ [meta.arm, ct], b, bi, (meta.platform ?: 'unknown') ]
            }
        }

        // Pool fragments of the same (sample_id, celltype) across lanes.
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

    // Determine if this run has per-cell BAMs (Smart-seq2) — affects celltype/goodCells strategy
    // Per-cell data uses SELECTIVE_MERGE; barcode-platform data uses sinto
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
    // Each BAM type is a labeled channel: [meta, bam, bai, bam_label]
    //
    // Cell filter × Dedup matrix:
    //                    none        coordinate(GATK)   umi(umi_tools)
    //   sample:          sample      sampleRmDup              sampleUmiDup
    //   goodCells:     goodCells goodCellsRmDup   goodCellsUmiDup
    //   annotated:       annotated   annotatedRmDup     annotatedUmiDup
    //   celltype:        celltype    celltypeRmDup      celltypeUmiDup
    //   recal:           —           recal (GATK + BQSR, qual=params.qual)
    //   recalq20:        —           recalq20 (GATK + BQSR, qual=20)
    //   recaldefault:    —           GATK-only, recal BAM + GATK defaults
    //                                (implicit whenever gatk + recal/recalq20)
    //
    // --dedup is comma-separated; all requested strategies run in parallel.
    // A tier is built only if it is in --bam_types (or is --percell_bam_type),
    // so sample-level runs never touch goodCells/annotated/celltype.

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

    // goodCells construction: annotated ∪ QC-passing unannotated
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

    // Per-cell platform samples (Smart-seq2). Empty in celltype-as-sample mode.
    ch_percell_src = toBool(params.split_by_celltype) ? Channel.empty() : ch_percell_samples

    // ── 3b. Coordinate-deduped BAMs (GATK MarkDuplicates) ──────────────
    //   Barcode platforms: MarkDuplicates on the merged BAM, cell-aware via CR.
    //   Per-cell platforms: dedup EACH cell, then merge (SS2 reads carry no cell
    //   tag, so a merged-BAM dedup would collapse fragments across cells). This
    //   merged BAM is sampleRmDup and feeds recal/recalq20.
    // Both dedup strategies (3b coordinate, 3c UMI) run per chromosome and per
    // GEM well (@RG SM; each well deduplicated separately), then merge
    // (modules/dedup.nf). A sample
    // that pools GEM groups WITHOUT @RG (teratoma merged arms) instead takes a
    // prebuilt BAM from the samplesheet (`rmdup_bam` / `umidup_bam` = per-
    // teratoma dedup, then merge). The prebuilt paths are read as a side
    // channel so `meta` -- and every cached task hash -- is unchanged.
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
    // [meta, chrom, file(s)...] per chromosome -> [meta, [files], ...] per sample,
    // emitted as soon as all dedup_chroms pieces of that sample are done.
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
    //   Same scatter / GEM-group / prebuilt routing as 3b (see above).
    if ('umi' in requested_dedup_strategies) {
        ch_umi_route = route_prebuilt.call(ch_base_bams, 'umidup_bam')
        UMI_DEDUP_CHROM(ch_umi_route.scatter.combine(Channel.fromList(dedup_chroms)), dedup_main)
        UMI_DEDUP_GATHER(scatter_gather_key.call(UMI_DEDUP_CHROM.out.bam))
        ch_umi_dedup = UMI_DEDUP_GATHER.out.bam.mix(ch_umi_route.prebuilt)
    }

    // Sample-level deduped labels
    ch_sampleRmDup  = ('coordinate' in requested_dedup_strategies)
                ? ch_coord_dedup.map { meta, bam, bai -> [ meta, bam, bai, 'sampleRmDup' ] }
                : Channel.empty()
    ch_sampleUmiDup = ('umi' in requested_dedup_strategies)
                ? ch_umi_dedup.map { meta, bam, bai -> [ meta, bam, bai, 'sampleUmiDup' ] }
                : Channel.empty()

    // ── 3d. Recalibrated BAMs (coordinate dedup + SplitNCigar + BQSR) ──
    // main_chroms handle (also (re)used by the caller stage below): restricts
    // GATK read traversal to requested chromosomes so externally-aligned BAMs
    // with extra/renamed scaffolds don't abort SplitNCigarReads/BaseRecalibrator.
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
        // recalq20: same BQSR'd BAM but variant calling uses qual=20 instead of params.qual
        ch_recalq20 = GATK_BQSR.out.bam.map { meta, bam, bai ->
            [ meta + [qual_override: 20], bam, bai, 'recalq20' ]
        }
        // recaldefault: GATK-only supplementary run — recal BAM, GATK defaults.
        // (GATK gets one extra run vs the other callers to show default-setting
        // performance; only GATK consumes this, not the other callers.)
        ch_gatk_default = GATK_BQSR.out.bam.map { meta, bam, bai -> [ meta, bam, bai, 'recaldefault' ] }
    } else {
        ch_recal = Channel.empty()
        ch_recalq20 = Channel.empty()
        ch_gatk_default = Channel.empty()
    }

    // ── 3e. goodCells — annotated ∪ QC-passing unannotated ───────────
    //        Barcode platforms: sinto filterbarcodes
    //        Per-cell platforms (SS2): selective samtools merge of matching BAMs

    // Split raw BAMs into barcode-platform vs per-cell-platform.
    // In celltype-as-sample mode each input is ALREADY one cell type, so the
    // secondary goodCells/annotated/celltype splits (3e/3f/3g, all derived from
    // ch_raw_barcode / ch_percell_src) must be disabled — emptying both
    // neutralizes every sinto/selective-merge call in that block.
    ch_raw_barcode = toBool(params.split_by_celltype)
        ? Channel.empty()
        : ch_raw_bams.filter { meta, bam, bai, annot, called_bc, cell_ids -> cell_ids == null }

    // Feed a tier's process only when one of its labels is wanted; an empty
    // input channel means the process is wired but never launches a task.
    // Per-cell platform RmDup tiers merge the per-cell-deduped BAMs from 3b.
    def gate = { List labels, ch -> labels.any { wants.call(it) } ? ch : Channel.empty() }

    // >>> ONLY BUILD THE CELL-TYPE TIER IF A CELLTYPE bam_type WAS ASKED FOR. <<<
    // SINTO_EXTRACT_CELLTYPE used to run unconditionally off ch_raw_barcode, so a
    // SAMPLE-LEVEL run was forced to split by cell type and then throw the result
    // away. That is not merely wasted compute: celltype_split.nf hard-errors when
    // an annotation maps one barcode to >1 cell type, which is UNAVOIDABLE for a
    // merged multi-sample BAM because 10x barcodes are reused across samples. The
    // guard is right — a raw-CB split of a merged BAM misassigns reads — but it
    // killed the whole run over a tier that was not requested. Measured: the
    // merged teratoma arms collide on 8 barcodes (1 AAVS1_ctrl, 7 NTC_ctrl) of
    // 14,511, and run 1710208 was cancelled by exactly that error on
    // merged_AAVS1_ctrl barcode ACCGTTCAGATAGCAT.
    //
    // Cell-type analysis is done properly by the celltype-as-sample track
    // (--split_by_celltype, which splits each lane before merging), so a
    // sample-level run has no reason to touch it.
    def wants_celltype = ['celltype', 'celltypeRmDup', 'celltypeUmiDup'].any { wants.call(it) }
    ch_ct_source = wants_celltype ? ch_raw_barcode : Channel.empty()
    if (!wants_celltype) {
        log.info "No celltype bam_type requested — skipping the cell-type split tier."
    }

    // Barcode platform → sinto
    ch_for_called_sinto = ch_raw_barcode
        .map { meta, bam, bai, annot, called_bc, cell_ids -> [ meta.id, meta, bam, bai ] }
        .combine(ch_good_bc, by: 0)
        .map { id, meta, bam, bai, bc -> [ meta, bam, bai, bc, 'goodCells' ] }
    SINTO_GOOD_CELLS(gate.call(['goodCells'], ch_for_called_sinto))
    ch_goodCells_sinto = SINTO_GOOD_CELLS.out.bam.map { meta, bam, bai ->
        [ meta, bam, bai, 'goodCells' ]
    }

    // Per-cell platform → selective merge
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
    //   Barcode platforms: sinto-filter the cell-aware MarkDuplicates BAM.
    //   Per-cell platforms: selective merge of the per-cell-deduped BAMs.
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
    // Barcode platform → sinto
    ch_for_annot_sinto = ch_raw_barcode.map { meta, bam, bai, annot, called_bc, cell_ids ->
        [ meta, bam, bai, annot, 'annotated' ]
    }
    SINTO_ANNOTATED(gate.call(['annotated'], ch_for_annot_sinto))
    ch_annotated_sinto = SINTO_ANNOTATED.out.bam.map { meta, bam, bai ->
        [ meta, bam, bai, 'annotated' ]
    }

    // Per-cell platform → selective merge using annotation file
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

    // Barcode platform → sinto celltype extraction
    SINTO_EXTRACT_CELLTYPE(gate.call(['celltype'], ch_ct_source.map { meta, bam, bai, annot, called_bc, cell_ids ->
        [ meta, bam, bai, annot ]
    }))
    ch_celltype_sinto = flattenCelltype(SINTO_EXTRACT_CELLTYPE.out.celltype_bams, 'celltype')

    // Per-cell platform → selective merge by celltype
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
    // Join strandedness with labeled BAMs for callers that need it
    // ch_strandedness is [meta, strandedness_string] keyed by meta.id
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
        // GATK also gets the supplementary default-settings run on the recal BAM
        // (recaldefault); other callers do not. Mixed only into GATK's input.
        GATK_HAPLOTYPECALLER(ch_all_labeled_bams.mix(ch_gatk_default), genome_fasta, main_chroms)
    }
    // Scatter units for the per-chromosome callers (scAllele; JACUSA2 with --jacusa2_scatter).
    def sc_chroms = params.main_chroms
        ? file(params.main_chroms).readLines().collect { it.trim() }.findAll { it }
        : file("${params.genome_fasta}.fai", checkIfExists: true).readLines().collect { it.tokenize('\t')[0] }
    def ch_jacusa2_variants = Channel.empty()
    if ('jacusa2' in requested_callers) {
        if (params.jacusa2_scatter) {
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
        // One task per chromosome, gathered in chromosome order (modules/variant_calling.nf).
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
    // Collect all variant outputs and run filtering on each
    snp_file      = params.snp_file     ? file(params.snp_file) : file('NO_FILE')
    strand_annots = params.strand_annots ? file(params.strand_annots) : file('NO_FILE')
    alu_annots    = params.alu_annots    ? file(params.alu_annots) : file('NO_FILE')

    // Each entry carries a per-sample genetic-variant blacklist path as a 5th
    // element. Genotyping callers (gatk 1/1, bcftools 1/1) emit one and are
    // joined by (id, bam_label); reditools3, jacusa2 and scallele pass NO_FILE.
    // scAllele's PASS is not a genotype call: its GLM was trained on genetic
    // variants with editing excluded, and well-supported edits PASS.
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

    // ── Benchmark axis: per-sample SNP subtraction (off = dbSNP only,
    // on = dbSNP + caller's own blacklist). 'on' only for callers that
    // produced a blacklist; pileup-only callers stay 'off' (NA).
    def snp_modes = params.snp_subtract_modes.tokenize(',')
    ch_variants_snpmodes = ch_all_variants.flatMap { meta, caller, bam_label, variants, blacklist ->
        def rows = []
        if ('off' in snp_modes) rows << [ meta, caller, bam_label, variants, 'off', no_blacklist ]
        if ('on'  in snp_modes && blacklist.name != 'NO_FILE') rows << [ meta, caller, bam_label, variants, 'on', blacklist ]
        return rows
    }

    // Expand with strand modes for benchmarking
    def requested_strand_modes = params.strand_modes.tokenize(',')
    ch_strand_modes = Channel.fromList(requested_strand_modes)

    // Which strand modes are meaningful for which caller.
    //
    // 'hybrid' only applies to bcftools and gatk, which carry I16/SB read-strand
    // counts. For the others hybrid would fall back to annotation and produce an
    // IDENTICAL duplicate, so those combos are dropped to avoid redundant work.
    //
    // 'caller' keeps the caller's OWN strand call and only applies to callers
    // that actually emit one: reditools3 (Strand column, 0/1/2/*), jacusa2, and
    // scallele (read-cluster strand of its RNA_EDITING? tag; caller mode keeps
    // only the tagged rows, see bin/scallele_io.py). bcftools and gatk hardcode
    // df$strand = "+", so 'caller' mode would be a constant "+" for them.
    //
    // NB 'caller' is deliberately NOT gated on library strandedness. On an
    // unstranded library the caller returns no strand and the mode degenerates
    // to "everything plus", misassigning mismatch types — which is precisely
    // what a user of that tool on that library would get, and is the finding
    // that argues for annotation-based stranding as the default.
    def hybrid_callers        = ['bcftools', 'gatk']
    def caller_strand_callers = ['reditools3', 'jacusa2', 'scallele']
    // Matched per MODE, not per caller. The previous form
    // (`strand_mode == 'annotation' || caller in hybrid_callers`) let ANY
    // non-annotation mode through for bcftools/gatk, so adding 'caller' to
    // --strand_modes would silently have produced constant-"+" runs for them.
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

    // scAllele goes through its own copy of DOWNSTREAM_R (modules/downstream_r.nf);
    // the ch_ds_* channels below merge both.
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
    // Add canonical, method-independent columns to the published filtered/pre100
    // tables using the SAME bam_type BAM that produced each call. When enabled,
    // DOWNSTREAM_R skips publishing the raw filtered/pre100 (REQUANT publishes
    // the requantified copies to the same 04_filtered path).
    if (toBool(params.requant_editing)) {
        // BAM lookup keyed by (sample id, bam_label); include GATK's recaldefault
        ch_bam_lut = ch_all_labeled_bams.mix(ch_gatk_default)
            .map { m, bam, bai, label -> [ [m.id, label], bam, bai ] }

        // Join filtered + pre100 on (id, caller, bam_label, strand_mode, snp_sub)
        ch_fp = ch_ds_filtered
            .map { m, c, bl, sm, ss, f -> [ [m.id, c, bl, sm, ss], m, c, bl, sm, ss, f ] }
            .join(
                ch_ds_pre100.map { m, c, bl, sm, ss, p -> [ [m.id, c, bl, sm, ss], p ] }
            )
            .map { key, m, c, bl, sm, ss, f, p -> [ [m.id, bl], m, c, bl, sm, ss, f, p ] }

        // Attach matching BAM (by id, bam_label), then strandedness -> libtype (by id)
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

    // Collate all stepwise mismatch tables into one file
    ch_ds_stepwise
        .map { meta, caller, bam_label, strand_mode, tsv -> tsv }
        .collectFile(
            name: 'all_stepwise_counts.tsv',
            storeDir: "${params.outdir}/summary",
            keepHeader: true,
            skip: 1
        )

    // Collate all per-run metrics into master results table
    ch_all_metrics = ch_ds_metrics
        .map { meta, caller, bam_label, strand_mode, tsv -> tsv }
        .collectFile(
            name: 'all_run_metrics.tsv',
            storeDir: "${params.outdir}/summary",
            keepHeader: true,
            skip: 1
        )

    // Generate side-by-side comparison report
    COLLATE_RESULTS(ch_all_metrics)

    // ── REDIportal overlap summary (human) ──────────────────────────────────
    // For every discovered site table, count known-editing (REDIportal) vs novel
    // sites + fractions -> summary/rediportal_overlap.tsv for figure plotting.
    // Gated on --rediportal so non-human datasets skip it. REDIportal is loaded
    // once; all pre100 tables are scored in a single R session.
    if (params.rediportal) {
        ch_pre100_all = ch_ds_pre100
            .map { m, c, bl, sm, ss, p -> p }
            .collect()
        REDIPORTAL_OVERLAP(ch_pre100_all, file(params.rediportal))
    }

    // Collate barcode coverage metrics from all filter processes
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

    // ── Step 6: Per-cell analyses ───────────────────────────────────────
    // sinto per-cell splitting only works on barcode platforms (CB tags).
    // Gated by params.percell_analyses (default true). Set false for caller-
    // comparison-only runs (e.g. scAllele-only): SINTO_SPLIT_PERCELL(pbmc2) is
    // memory-heavy and irrelevant to site-level caller benchmarking.
    if (toBool(params.percell_analyses)) {
    ch_barcode_base_bams = ch_base_bams
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .combine(
            ch_barcode_samples.map { meta, bams, bais, annot, called_bc, cell_ids -> [ meta.id ] },
            by: 0
        )
        .map { id, meta, bam, bai -> [ meta, bam, bai ] }

    SINTO_SPLIT_PERCELL(ch_barcode_base_bams.map { m, b, bi -> [ m, b, bi, file('NO_FILE') ] })
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
        // Use the BAM channel matching the requested BAM type
        // ch_base_bams is the merged BAM (sample equivalent)
        // For other BAM types, route to the appropriate channel
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

        // Per-sample barcode list restricting which cells to split (avoids
        // splitting background barcodes when using the sample BAM). In
        // split_by_celltype mode each BAM is one cell type already, so all of
        // its CBs are used (NO_FILE).
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

        // Split the chosen BAM into per-cell BAMs ONCE (sorted+indexed), then
        // discover per cell -- avoids the O(cells x whole-BAM) re-scan.
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
