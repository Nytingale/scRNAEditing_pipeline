# scRNAEditing

Nextflow pipeline for benchmarking A-to-I RNA editing detection in single-cell RNA-seq transcriptomes

## Overview

This pipeline systematically runs variant callers (BCFtools, REDItools3, GATK HaplotypeCaller, JACUSA2 and scAllele) across multiple BAM processing strategies for detecting adenosine-to-inosine (A-to-I) RNA editing in single-cell transcriptomes. Users may run any combination of the method configurations available.

### Pipeline steps

1. **Input validation** — Check BAM tags (CR, CB, UB, MD, NH), add MD tag if missing, filter for uniquely-mapped reads
2. **Optional remapping** — STARsolo alignment with platform-aware parameters for different sequencing chemistries (10x V2, 10x V3, Smart-seq2, DroNc-seq, Drop-seq, Cel-seq, sci-seq, InDrops)
3. **Strandedness inference** — RSeQC infer_experiment.py
4. **BAM type generation** — sample / goodCells / annotated / celltype × unprocessed / RmDup / UmiDup, plus recal / recalq20 / recalDefault
5. **Variant calling** — All requested caller x BAM type combinations
6. **Downstream filtering** — Caller's own strand (only for REDItools3 and JACUSA2), Gene annotation-guided strand assignment otherwise, SNP removal with known list (or also with per sample blacklist only available if using variant callers BCFtools, GATK or scAllele), adaptive window optimization for enrichment of A-to-G mismatches
7. **Optional per-cell analysis** — Per-cell variant calling and Alu Editing Index (AEI)

## Quick start

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/genome.fa \
    --genome_gtf /path/to/genes.gtf \
    --snp_file /path/to/snp_GRanges.RData \
    --strand_annots /path/to/strand_annotations.RData \
    --main_chroms /path/to/main_chroms.txt \
    --callers bcftools,reditools3 \
    --bam_types sample,sampleRmDup \
    --dedup coordinate \
    -profile conda
```

## Samplesheet format

CSV with header. See `assets/samplesheet_example.csv`.

| Column | Required | Description |
|---|---|---|
| `sample_id` | Yes | Unique sample identifier |
| `bam` | Yes* | Path to mapped BAM file |
| `celltype_annot` | Yes | 2-column TSV (no header): barcode, celltype |
| `called_barcodes` | No | Fallback cell list when no count matrix is available (Smart-seq2 QC list) |
| `count_matrix` | No | STARsolo/CellRanger matrix directory; auto-discovered from `<bam>_Solo.out` if omitted |
| `qc_barcodes` | No | Precomputed goodCells barcode list (skips metric QC) |
| `cell_id` | Smart-seq2 | Per-cell BAM identifier used by selective merge |
| `platform` | If remapping | One of: 10xv2, 10xv3, dropseq, celseq2, seqwell, dronc, sciseq, indrop, indropv2, smartseq2 |
| `fastq_r1` | If remapping | R1 FASTQ (barcode read) |
| `fastq_r2` | If remapping | R2 FASTQ (cDNA read) |

*Provide either `bam` (BAM input mode) or `fastq_r1` + `fastq_r2` + `platform` (remap mode).

## Parameters

### Required

| Parameter | Description |
|---|---|
| `--input` | Samplesheet CSV |
| `--genome_fasta` | Reference genome FASTA |
| `--genome_gtf` | Gene annotation GTF |

### Genome FASTA notes

The provided annotation files (gene/Alu GRanges, SNP positions) are **coordinate-based** and work with any hg38/GRCh38 variant — whether labelled "hg38", "GRCh38", or "GRCh38_analysis_set" — as long as chromosomes use UCSC naming (`chr1`, `chr2`, ..., `chrX`, `chrY`). No separate annotation build is needed.

**Ambiguous bases:** Some genome FASTA files contain IUPAC ambiguity codes (R, B, Y, etc.) at certain positions. These can cause false positive mismatches with some variant callers (notably older versions of REDItools). We provide an optional utility to convert non-ACTG bases to N:

```bash
bash nextflow/bin/mask_non_actg.sh input.fa > masked.fa
samtools faidx masked.fa
```

This masking is optional and may not be needed with newer tool versions. If in doubt, it is safer to mask.

### Recommended

| Parameter | Default | Description |
|---|---|---|
| `--snp_file` | null | Known SNP positions — RData with `snp_gr` GRanges (preferred) or BED file |
| `--known_sites_vcf` | null | Known sites VCF/BED for GATK BQSR recalibration (separate from SNP filtering) |
| `--strand_annots` | null | RData file with `plus.annot` and `minus.annot` GRanges objects |
| `--alu_annots` | null | RData file with `alu.plus.annot` and `alu.minus.annot` GRanges (human only) |
| `--main_chroms` | null | One chromosome name per line; restricts analysis to these chromosomes |
| `--refseq_bed` | null | RefSeq BED12 for strandedness inference |
| `--repeat_annots` | null | Repeat element BED for per-cell AEI calculation |

### Processing options

| Parameter | Default | Description |
|---|---|---|
| `--remap` | false | Re-align from FASTQs with STARsolo |
| `--dedup` | coordinate | Deduplication: `none`, `coordinate` (GATK MarkDuplicates with CR tag), or `umi` (umi_tools) |
| `--callers` | bcftools,reditools3,gatk,jacusa2 | Comma-separated list of variant callers (`scallele` is opt-in: `--callers ...,scallele`) |
| `--scallele_run_mode` | Variant_Caller | scAllele `--run_mode`: `Variant_Caller`, `Full`, or `Training` |
| `--scallele_conda` | `/home/yng7/miniforge3/envs/scAllele` | Conda env for SCALLELE only (do not mix into scRNAEditing) |
| `--snp_subtract_modes` | `off,on` | Benchmark axis for per-sample SNP filtering. `off` = dbSNP only; `on` = dbSNP + the caller's own per-sample genetic-variant blacklist. Only genotyping-capable callers (scallele PASS, gatk/bcftools homozygous-ALT `1/1`) run `on`; reditools3/jacusa2 are `off` only (NA). Outputs/metrics carry a `snp_subtract` label; publishDir gains a `snp_{off,on}/` level. |
| `--bam_types` | sample,sampleRmDup,sampleUmiDup,<br>recal,recalq20,<br>goodCells,goodCellsRmDup,goodCellsUmiDup,<br>annotated,annotatedRmDup,annotatedUmiDup,<br>celltype,celltypeRmDup,celltypeUmiDup | Comma-separated list of BAM processing types |
| `--qual` | 30 | Minimum base/mapping quality |
| `--min_cov` | 10 | Minimum read coverage at a site |
| `--min_alt` | 3 | Minimum edited (alt) read count at a site, under `--alt_rule count`. See [Alt-read rule](#alt-read-rule-and-the-sensitivity-floor) — this, not `--min_edit_level`, is the real sensitivity threshold |
| `--alt_rule` | count | How a site qualifies as edited: `count` (alt reads >= `--min_alt`) or `binom` (alt count significantly above `--seq_error_rate`) |
| `--seq_error_rate` | 0.001 | Per-base error under the null, `--alt_rule binom` only |
| `--alt_binom_alpha` | 0.05 | Significance threshold, `--alt_rule binom` only |
| `--min_alt_floor` | 2 | Never accept a single alt read, `--alt_rule binom` only |
| `--adaptive_windows` | 0,20,40,...,200 | Comma-separated window sizes for adaptive filtering |
| `--star_index` | null | Pre-built STAR genome index directory (built automatically if not provided) |
| `--outdir` | results | Output directory |
| `--mito_chr` | auto from GTF | Mitochondrial chromosome name(s), comma-separated; override auto-detect |
| `--qc_min_features` | 200 | goodCells: min genes detected in unannotated cells |
| `--qc_max_features` | 6000 | goodCells: max genes detected in unannotated cells |
| `--qc_max_pct_mt` | 20 | goodCells: max mitochondrial percent |
| `--qc_min_counts` | 500 | goodCells: min UMIs/counts in unannotated cells |

### Alt-read rule and the sensitivity floor

`--min_alt`, not `--min_edit_level`, is what actually limits sensitivity. Of the
three filters in the coverage/editing step, only one does any work at typical
10x depth:

| threshold | how often it binds |
|---|---|
| `--min_edit_level 0.01` | 2.6% of sites (needs >300 reads to bind) — effectively inert |
| `--min_cov 10` | excludes ~14% — mild |
| **`--min_alt 3`** | **binds essentially always** |

Because `--min_alt 3` implies `editing >= 3/coverage`, the effective editing
floor is set by depth: ~19% at the 25th coverage percentile, ~8% at the median,
~2% at the 90th. At well-covered sites (>=50 reads) the median editing level is
only ~5%, so at median depth the 3-read rule discards most of the low-level
editing that is visible where depth allows.

**If you care about low-level editing, lower `--min_alt` to 2.** That is the
supported way to do it.

> **`--alt_rule binom` does not do what it looks like it does.** It was added to
> replace the fixed count with a test against sequencing error, on the
> expectation that it would recover low-level sites while still rejecting
> single-read noise. Measured on 2.8M real reditools3 rows, it is **exactly
> equivalent to `count` with `--min_alt 2`** — identical surviving site counts,
> in both strand modes. With `--min_cov 10` and `--min_alt_floor 2`, two edited
> reads are already significant against a 0.1% error rate at any coverage the
> pipeline admits (p ~ 4e-5 at coverage 10), so the floor decides every case and
> the test never binds.
>
> Correcting for multiple testing does not rescue it: the rule evaluates ~2.4M
> sites, and at a Bonferroni-corrected alpha it requires ~5 edited reads at
> median depth — *stricter* than the `--min_alt 3` it was meant to relax.
>
> `binom` is kept because it becomes meaningful at much higher per-site depth
> than 10x scRNA-seq provides, where the test can bind before the floor does.
> **For typical single-cell data, prefer `--alt_rule count --min_alt 2`**, which
> is clearer about what it is doing. `count` remains the default.

## BAM type hierarchy

Cell filter × dedup suffix. `goodCells` is built **before** BAM filtering:

1. Candidate pool = STARsolo/CellRanger filtered barcodes ∪ annotated barcodes
2. `unannotated = pool − annotated` (annotated cells are assumed QC-passed)
3. QC unannotated cells only: nFeature, nCount, pct_mt (mito genes from GTF chromosome)
4. `goodCells = annotated ∪ QC-passing unannotated`

| Transition | Isolated effect |
|---|---|
| `sample` → `goodCells` | cell quality (empty droplets / low-quality removed) |
| `goodCells` → `annotated` | annotation completeness |
| `annotated` → `celltype` | cell-type resolution |

Without a count matrix (Smart-seq2), the samplesheet `cell_id` list / `called_barcodes` is treated as already QC-filtered. |

## Deduplication strategies

| Strategy | Flag | Tool | Use case |
|---|---|---|---|
| Coordinate-based | `--dedup coordinate` | GATK MarkDuplicates (--BARCODE_TAG CR) | Default; works for all protocols |
| UMI-aware | `--dedup umi` | umi_tools dedup (--per-cell --cell-tag=CB --umi-tag=UB) | 10x, Drop-seq, and other UMI protocols |
| None | `--dedup none` | Passthrough | Skip deduplication entirely |

## Supported platforms (for --remap)

| Platform | CB length | UMI length | Whitelist required |
|---|---|---|---|
| 10xv2 | 16 | 10 | 737K-august-2016.txt |
| 10xv3 | 16 | 12 | 3M-february-2018.txt |
| dropseq | 12 | 8 | No |
| celseq2 | 6 | 6 | No |
| seqwell | 12 | 8 | No |
| dronc | 12 | 8 | No |
| sciseq | 10 | 8 | No |
| indrop | 8 | 6 | No |
| indropv2 | Complex | Complex | User-supplied |
| smartseq2 | N/A | N/A | N/A (plate-based) |

Place whitelist files in `assets/whitelists/`.

## Output structure

```
results/
├── {sample_id}/
│   ├── 00_validated/           # Validated, MD-tagged, uniquely mapped BAMs
│   ├── 00_qc/                  # Strandedness reports
│   ├── 01_starsolo/            # STARsolo output (if --remap)
│   ├── 02_dedup/               # Deduplicated BAMs + metrics
│   ├── 02_goodcells/           # goodCells barcode list + QC summary
│   ├── 02_barcode_filtered/    # goodCells / annotated filtered BAMs
│   ├── 02_gatk_preprocess/     # SplitNCigar + BQSR BAMs
│   ├── 03_variants/
│   │   ├── bcftools/{bam_type}/    # bcftools mpileup output
│   │   ├── reditools3/{bam_type}/  # REDItools3 output
│   │   ├── gatk/{bam_type}/        # GATK HaplotypeCaller output
│   │   ├── jacusa2/{bam_type}/     # JACUSA2 call-1 output
│   │   └── scallele/{bam_type}/    # scAllele SNV TSV (opt-in caller)
│   ├── 03_celltype_bams/       # Per-celltype BAMs
│   ├── 04_filtered/
│   │   └── {caller}/{bam_type}/    # Filtered sites + summary stats
│   ├── 04_percell_bams/        # Per-cell BAMs
│   ├── 05_percell_variants/    # Per-cell variant calls
│   └── 05_percell_aei/         # Per-cell Alu Editing Index
├── pipeline_trace.txt          # Task-level resource usage
├── pipeline_timeline.html      # Interactive timeline
├── pipeline_report.html        # Execution report
└── pipeline_dag.html           # Workflow DAG
```

## Execution profiles

```bash
-profile standard    # Local execution
-profile sge         # Sun Grid Engine
-profile pbs         # PBS/Torque
-profile slurm       # SLURM
-profile conda       # Conda environment
-profile docker      # Docker container
-profile singularity # Singularity container
```

Profiles can be combined: `-profile slurm,conda`

## Resource tracing

The pipeline automatically records CPU, memory, and wall time for every task via Nextflow's built-in tracing. The `pipeline_trace.txt` file contains per-task metrics including peak RSS, peak virtual memory, read/write bytes, and duration. This directly addresses reproducibility and resource documentation requirements.

## Software requirements

- Nextflow >= 23.04
- samtools
- bcftools
- STAR (if remapping)
- GATK 4
- JACUSA2
- REDItools3 (pip install reditools)
- umi_tools (if --dedup umi)
- sinto
- RSeQC
- R >= 4.0 with packages: optparse, data.table, parallel, GenomicRanges, janitor, Rcpp, stringr

## Differential editing (separate)

After inspecting the discovered editing sites, differential editing analysis can be performed using the standalone scripts in `../diff_editing/`. See `diff_editing/README.md` for usage. This includes REDIT-LLR, Wilcoxon, and limma-based tests with a union consensus approach.
