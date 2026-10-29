# scRNAEditing

> [!WARNING]
> **This pipeline is still in development and is not ready for use.** This repository accompanies a manuscript under revision. The pipeline, analysis scripts and documentation are still being tested and changed, and parts may not yet run outside our environment. Please wait for a tagged release before downloading it for your own analyses.

Nextflow pipeline for benchmarking A-to-I RNA editing detection in single-cell RNA-seq transcriptomes

## Overview

This pipeline runs variant callers (BCFtools, REDItools3, GATK HaplotypeCaller, JACUSA2, scAllele and Cellsnp-lite) across multiple BAM processing strategies to detect adenosine-to-inosine (A-to-I) RNA editing, seen as A-to-G mismatches, in single-cell transcriptomes. Users may run any combination of the method configurations available. A method configuration is one end-to-end combination of caller, BAM type, strand mode and per-sample SNP removal.

### Pipeline steps

1. **Input validation**: check BAM tags (CR, CB, UB, MD, NH), add the MD tag if missing, keep uniquely mapped reads
2. **Optional remapping**: STARsolo alignment with platform-aware parameters (10x v2, 10x v3, Smart-seq2, DroNc-seq, Drop-seq, CEL-seq2, Seq-Well, sci-seq, inDrops)
3. **Lane merging**: rows sharing a `sample_id` are lanes of one sample and are merged before any analysis
4. **Strandedness inference**: RSeQC `infer_experiment.py`
5. **BAM type generation**: sample / goodCells / annotated / celltype × no dedup / coordinate dedup (RmDup) / UMI dedup (UmiDup), plus recal / recalq20 (BQSR) and recaldefault (GATK only)
6. **Variant calling**: every requested caller × BAM type. Deduplication, scAllele and (optionally) JACUSA2 run as one task per chromosome, then gather, so deep BAMs stay within memory
7. **Downstream filtering**: parsing, main chromosomes, removal of multi-mismatch sites, known-SNP removal (optionally with a per-sample genetic-variant blacklist, BCFtools / GATK only), strand assignment, coverage and editing-level filters, and optional removal of sites near other mismatch types
8. **Re-quantification**: one pileup-based, strand-aware editing level and coverage for every published site, comparable across callers
9. **Summaries**: per-run metrics, stepwise mismatch counts and, for human data, REDIportal overlap
10. **Optional per-cell analysis**: per-cell variant calling and Alu Editing Index (AEI)

## Installation

All software the pipeline runs is in one container image, `ghcr.io/nytingale/scrnaediting:1.0.0`. It holds the four conda environments the pipeline uses, rebuilt from the exact lock files in `envs/`:

| Environment | Contents | Why separate |
|---|---|---|
| `scRNAEditing` | Nextflow, samtools, BCFtools, GATK, JACUSA2, REDItools3, STAR, umi_tools, sinto, RSeQC, R and all R packages | main environment |
| `scAllele` | scAllele 0.0.9.7 | needs numpy 1.26 / scipy 1.15; the main environment has numpy 2 |
| `cellsnp` | Cellsnp-lite 1.2.3 | built against htslib 1.24; the main environment's samtools uses htslib 1.23.1 |
| `cutadapt` | cutadapt 4.9, pigz | Smart-seq2 TSO trimming (`--tso_clip`, `--remap`); needs Python 3.12, the main environment has 3.11 |

Use the container (Singularity / Apptainer on HPC clusters, or Docker) unless you cannot; the conda route below reproduces the same environments on linux-64.

### 1. Nextflow and the pipeline

Nextflow needs Java 17 or later. The benchmark used Nextflow 26.04.6:

```bash
curl -s https://get.nextflow.io | bash        # installs ./nextflow
export NXF_VER=26.04.6                        # optional: pin the version used for the benchmark
git clone https://github.com/Nytingale/scRNAEditing_pipeline.git
cd scRNAEditing_pipeline
```

### 2a. Singularity or Apptainer (recommended on clusters)

Add `-profile singularity` (or `apptainer`) to every run, together with your scheduler profile:

```bash
nextflow run main.nf -profile singularity,slurm --input samplesheet.csv ...
```

Nextflow downloads the image and converts it once to a `.img` file. Set a shared cache so this happens only once per cluster; the cache directory must be on a file system the compute nodes can read:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/containers     # NXF_APPTAINER_CACHEDIR for Apptainer
```

If compute nodes have no internet access, pull the image on a login node into that directory first, under the file name Nextflow expects:

```bash
singularity pull "$NXF_SINGULARITY_CACHEDIR/ghcr.io-nytingale-scrnaediting-1.0.0.img" \
    docker://ghcr.io/nytingale/scrnaediting:1.0.0
```

or build/pull a `.sif` anywhere and pass it directly with `--container /path/to/scrnaediting_1.0.0.sif`. Check the image with:

```bash
singularity exec "$NXF_SINGULARITY_CACHEDIR/ghcr.io-nytingale-scrnaediting-1.0.0.img" \
    bash -c 'samtools --version | head -1; bcftools --version | head -1; scAllele --help | head -1; cellsnp-lite --version'
```

`singularity.autoMounts` is on, so input paths outside your home directory are bound automatically. If your cluster needs extra bind paths, add `singularity.runOptions = '-B /data'` in a `-c` config file.

### 2b. Docker

```bash
docker pull ghcr.io/nytingale/scrnaediting:1.0.0
nextflow run main.nf -profile docker --input samplesheet.csv ...
```

Tasks run as your user (`-u $(id -u):$(id -g)`), so outputs are not owned by root.

### 2c. Conda (no containers)

Create the environments from the lock files. The lock files list every package by exact build and URL, so no solving takes place; they are for linux-64 only. `envs/*.yml` list the same versions without build strings for other platforms, but are not guaranteed to solve identically.

```bash
# main environment, then the R packages that were installed from CRAN / GitHub after conda
mamba create -n scRNAEditing --file envs/scRNAEditing.linux-64.lock.txt
mamba run -n scRNAEditing Rscript envs/install_r_packages.R

# scAllele: conda packages, then the exact PyPI packages
mamba create -n scAllele --file envs/scAllele.linux-64.lock.txt
mamba run -n scAllele pip install --no-deps -r envs/scAllele.pip.txt

# Cellsnp-lite
mamba create -n cellsnp --file envs/cellsnp.linux-64.lock.txt

# cutadapt (Smart-seq2 TSO trimming with --remap)
mamba create -n cutadapt --file envs/cutadapt.linux-64.lock.txt
```

`install_r_packages.R` installs 15 packages at the versions used (`envs/scRNAEditing.r_packages.tsv`, e.g. Seurat 5.5.1, data.table 1.18.6.1, DoubletFinder from a fixed GitHub commit) and stops if any version differs. It compiles from source with the compilers in the environment, which took about 6 minutes in our test.

Then point the pipeline at the environments:

```bash
CONDA_BASE=$(conda info --base)
nextflow run main.nf -profile conda,slurm \
    --conda_env      "$CONDA_BASE/envs/scRNAEditing" \
    --scallele_conda "$CONDA_BASE/envs/scAllele" \
    --cellsnp_bin    "$CONDA_BASE/envs/cellsnp/bin/cellsnp-lite" \
    --cutadapt_conda "$CONDA_BASE/envs/cutadapt" \
    --input samplesheet.csv ...
```

Instead of `--conda_env` you can `export SCRNA_CONDA_ENV=$CONDA_BASE/envs/scRNAEditing`. Per-cell scAllele (`--percell_callers scallele`) activates the environment by name, so keep the name `scAllele`.

Without `--conda_env`, Nextflow builds the main environment itself from the lock file (in `work/conda/`). That environment has the conda versions of the 15 R packages above, not the versions used for the benchmark, so prefer the explicit environments for exact reproduction.

### Building the image yourself

`containers/Dockerfile` builds the image from `envs/`; `.github/workflows/container.yml` builds and publishes it for every release tag. To build locally (Docker or Podman), from the repository root:

```bash
docker build -f containers/Dockerfile -t scrnaediting:1.0.0 .
singularity build scrnaediting_1.0.0.sif docker-daemon://scrnaediting:1.0.0
```

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
    --dedup none,coordinate \
    -profile singularity
```

## Recommended configuration

In the teratoma ADAR1 knockout benchmark (pseudobulk; ground truth from sites with decreased editing in the knockout), the top-ranked configuration was BCFtools on deduplicated BAMs (coordinate or UMI dedup), `hybrid` strand mode and the per-sample SNP blacklist (`snp_on`), with removal of sites near other mismatch types:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/genome.fa \
    --genome_gtf /path/to/genes.gtf \
    --snp_file /path/to/snp_GRanges.RData \
    --strand_annots /path/to/strand_annotations.RData \
    --main_chroms /path/to/main_chroms.txt \
    --callers bcftools \
    --bam_types sampleRmDup,sampleUmiDup \
    --dedup coordinate,umi \
    --strand_modes hybrid \
    --snp_subtract_modes on \
    --mismatch_window_filter true \
    -profile singularity
```

Best configuration of each other caller, ranked across the knockout truth sets: JACUSA2 5th, REDItools3 7th, scAllele 9th, GATK 27th, Cellsnp-lite 32nd.

## Reference files

The reference files used for the human benchmark (GRCh38 analysis set):

| Parameter | File | Origin |
|---|---|---|
| `--genome_fasta` | `hg38.analysisSet.fa` | UCSC hg38 analysis set, no ALT contigs |
| `--genome_gtf` | `hg38.ncbiRefSeq.gtf` | UCSC NCBI RefSeq track |
| `--main_chroms` | `chrlist_main_chrs.txt` | chr1-22, X, Y |
| `--snp_file` | `dbSNP156_common_GRanges.RData` | **built**: common SNPs of dbSNP 156 as GRanges |
| `--strand_annots`, `--alu_annots` | `extended_annots_GRCh38_analysis_set.RData` | **built**: NCBI RefSeq annotation of the GRCh38 analysis set, per strand, extended with introns and intergenic regions, plus UCSC RepeatMasker Alu elements |
| `--known_sites_vcf` | `dbSNP156_GRCh38.vcf.gz` | dbSNP 156, chromosome names converted to UCSC |
| `--rediportal` | `REDIportal_hg38_v2.vcf.gz` | REDIportal v2 converted to VCF |
| `--refseq_bed` | `hg38_RefSeq.bed` | RSeQC BED12 |
| `--star_index` | `stardb_hg38_<overhang>/` | STAR `genomeGenerate` (built by the pipeline if not given) |

Get them in one of two ways:

1. **Download the built files** (extended annotations and SNP GRanges for GRCh38, mm39 and six other genomes, and the REDIportal VCFs for human and mouse) from Zenodo with `bash references/download_reference_files.sh <dir>`, and the public files with the download scripts below.
2. **Rebuild everything** with the scripts in `references/` (`human/01` to `08`, `mouse/01` to `07`; they write to `$REFDIR`). The rebuilt files match the ones used for the paper, except the human gene annotation: NCBI has since replaced its GTF, so a rebuild differs in about 0.5% of rows. Use the Zenodo file to reproduce the paper exactly.

`references/README.md` describes each step, the sources and versions, the other genomes, and how to build the annotations for a new genome.

## Samplesheet format

CSV with header. See `assets/samplesheet_example.csv` and `SAMPLESHEET_SPEC.md`.

| Column | Required | Description |
|---|---|---|
| `sample_id` | Yes | Biological sample. Rows sharing a `sample_id` are lanes of that sample and are merged |
| `bam` | Yes* | Path to mapped BAM file |
| `celltype_annot` | Yes | 2-column TSV (no header): barcode, celltype |
| `called_barcodes` | No | Fallback cell list when no count matrix is available (Smart-seq2 QC list) |
| `count_matrix` | No | STARsolo / Cell Ranger matrix directory; auto-discovered from `<bam>_Solo.out` if omitted |
| `qc_barcodes` | No | Precomputed goodCells barcode list (skips metric QC) |
| `rmdup_bam`, `umidup_bam` | No | Prebuilt coordinate / UMI-deduplicated BAM for samples that pool GEM wells without `@RG` tags |
| `cell_id` | Smart-seq2 | Per-cell BAM identifier used by selective merge |
| `platform` | If remapping | One of: 10xv2, 10xv3, dropseq, celseq2, seqwell, dronc, sciseq, indrop, indropv2, smartseq2 |
| `fastq_r1` | If remapping | R1 FASTQ (barcode read) |
| `fastq_r2` | If remapping | R2 FASTQ (cDNA read) |

*Provide either `bam` (BAM input mode) or `fastq_r1` + `fastq_r2` + `platform` (remap mode).

Distinct biological samples need distinct `sample_id`s: 10x barcodes are reused across samples, so a merged multi-sample BAM cannot be split by barcode. To pool cell types across samples, use `--split_by_celltype`, which splits each lane by cell type before merging.

## Parameters

### Required

| Parameter | Description |
|---|---|
| `--input` | Samplesheet CSV |
| `--genome_fasta` | Reference genome FASTA (indexed: `.fai`) |
| `--genome_gtf` | Gene annotation GTF |

### Software

| Parameter | Default | Description |
|---|---|---|
| `--container` | `ghcr.io/nytingale/scrnaediting:1.0.0` | Container profiles: image tag, digest or local `.sif` |
| `--conda_env` | `$SCRNA_CONDA_ENV`, else `envs/scRNAEditing.linux-64.lock.txt` | conda profile: main environment (prefix or lock file) |

### Genome FASTA notes

The provided annotation files (gene/Alu GRanges, SNP positions) are **coordinate-based** and work with any hg38/GRCh38 variant — whether labelled "hg38", "GRCh38", or "GRCh38_analysis_set" — as long as chromosomes use UCSC naming (`chr1`, `chr2`, ..., `chrX`, `chrY`). No separate annotation build is needed.

**Ambiguous bases:** Some genome FASTA files contain IUPAC ambiguity codes (R, B, Y, etc.) at certain positions. These can cause false positive mismatches with some variant callers (notably older versions of REDItools). We provide an optional utility to convert non-ACTG bases to N:

```bash
bash bin/mask_non_actg.sh input.fa > masked.fa
samtools faidx masked.fa
```

This masking is optional and may not be needed with newer tool versions. If in doubt, it is safer to mask. Note that this is NOT soft or hard masking repetitive regions.

### Recommended

| Parameter | Default | Description |
|---|---|---|
| `--snp_file` | null | Known SNP positions: RData with `snp_gr` GRanges (preferred) or BED file |
| `--known_sites_vcf` | null | Known sites VCF/BED for GATK BQSR (separate from SNP filtering; null = no recalibration) |
| `--strand_annots` | null | RData file with `plus.annot` and `minus.annot` GRanges objects |
| `--alu_annots` | null | RData file with `alu.plus.annot` and `alu.minus.annot` GRanges (human only) |
| `--main_chroms` | null | One chromosome name per line; restricts analysis to these chromosomes |
| `--refseq_bed` | null | RefSeq BED12 for strandedness inference |
| `--rediportal` | null | REDIportal VCF (human); writes `summary/rediportal_overlap.tsv` (known vs novel sites) |
| `--repeat_annots` | null | Repeat element BED for per-cell AEI calculation |

### BAM processing

| Parameter | Default | Description |
|---|---|---|
| `--remap` | false | Re-align from FASTQs with STARsolo |
| `--tso_clip` | true | With `--remap`, remove the template-switch oligo (TSO) before alignment. 10xv2 / 10xv3: STAR `--clipAdapterType CellRanger4` (5' TSO and 3' poly(A), as Cell Ranger >= 4). smartseq2: cutadapt trims the TSO from the 5' end of both mates (`CUTADAPT_TSO`; STAR's CellRanger4 mode clips mate 1 only). Without it, the TSO's final GGG aligns onto genomic G-A-G and creates A-to-G artifacts |
| `--tso_seq` | `AAGCAGTGGTATCAACGCAGAGTACATGGG` | TSO sequence for the Smart-seq2 trim (10x and Smart-seq2 share it) |
| `--tso_min_overlap` | 5 | Shortest TSO 3' end (`ATGGG`) trimmed at a read start (Smart-seq2) |
| `--tso_min_length` | 20 | Smart-seq2: drop read pairs with a mate shorter than this after trimming |
| `--cutadapt_conda` | `envs/cutadapt.yml` | conda profile only: cutadapt environment (prefix or file) for `CUTADAPT_TSO`; ignored with containers |
| `--star_index` | null | Pre-built STAR genome index directory (built automatically if not provided) |
| `--bam_types` | sample,sampleRmDup,sampleUmiDup,<br>recal,recalq20,<br>goodCells,goodCellsRmDup,goodCellsUmiDup,<br>annotated,annotatedRmDup,annotatedUmiDup,<br>celltype,celltypeRmDup,celltypeUmiDup | Comma-separated list of BAM types; only requested types are built |
| `--dedup` | `none,coordinate` | Comma-separated dedup strategies: `none`, `coordinate` (GATK MarkDuplicates), `umi` (umi_tools). RmDup types need `coordinate`, UmiDup types need `umi` |
| `--dedup_cell_tag` | `CR` | Cell barcode tag for dedup grouping: `CR` (raw; recommended) or `CB` (corrected; reads with `CB:Z:-` are dropped) |
| `--split_by_celltype` | false | Each cell type becomes its own sample (`<sample>__<celltype>`) and runs the full sample-level matrix; use sample-level `--bam_types` |
| `--min_celltype_cells` | 0 | With `--split_by_celltype`, drop cell types with fewer cells (0 = keep all) |
| `--mito_chr` | auto from GTF | Mitochondrial chromosome name(s), comma-separated; overrides auto-detection |
| `--qc_min_features` | 200 | goodCells: min genes detected in unannotated cells |
| `--qc_max_features` | 6000 | goodCells: max genes detected in unannotated cells |
| `--qc_max_pct_mt` | 20 | goodCells: max mitochondrial percent |
| `--qc_min_counts` | 500 | goodCells: min UMIs/counts in unannotated cells |
| `--outdir` | results | Output directory |

### Variant calling

| Parameter | Default | Description |
|---|---|---|
| `--callers` | bcftools,reditools3,gatk,jacusa2 | Comma-separated list of callers; `scallele` and `cellsnp` are opt-in |
| `--qual` | 30 | Minimum base and mapping quality (recalq20 uses 20) |
| `--bcftools_baq` | true | bcftools mpileup base alignment quality (BAQ). `true` keeps the bcftools default (partial BAQ); `false` passes `-B` (BAQ off). In the teratoma ground truth, BAQ off lost 0.02-0.10 F1 on deep (merged) data and gained <= 0.01 on shallow data |
| `--jacusa2_min_score` | 1.25 | Minimum JACUSA2 score |
| `--jacusa2_scatter` | false | Run JACUSA2 as one task per chromosome (for very deep BAMs) instead of one whole-genome task |
| `--scallele_conda` | `envs/scAllele.yml` | conda profile only: scAllele environment (prefix or file); ignored with containers |
| `--scallele_run_mode` | Variant_Caller | scAllele `--run_mode`: `Variant_Caller`, `Full`, or `Training` |
| `--cellsnp_bin` | `cellsnp-lite` | Cellsnp-lite binary; with the conda profile, the one in the `cellsnp` environment |

### Downstream filtering

| Parameter | Default | Description |
|---|---|---|
| `--strand_modes` | annotation | Comma-separated strand modes to benchmark; see [Strand modes](#strand-modes) |
| `--strand_ambiguity_threshold` | 0.3 | `hybrid` mode: read-strand ambiguity cutoff |
| `--snp_subtract_modes` | `off,on` | Per-sample SNP removal. `off` = known SNPs only; `on` = known SNPs + the caller's own per-sample genetic-variant blacklist (homozygous-ALT calls). Only BCFtools and GATK run `on`; other callers are `off` only. Outputs gain a `snp_{off,on}/` level |
| `--min_cov` | 10 | Minimum read coverage at a site (all callers) |
| `--min_alt` | 3 | Minimum edited (alt) read count, under `--alt_rule count`. See [Alt-read rule](#alt-read-rule-and-the-sensitivity-floor): this, not `--min_edit_level`, is the real sensitivity threshold |
| `--min_edit_level` | 0.01 | Minimum editing level |
| `--alt_rule` | count | How a site qualifies as edited: `count` (alt reads >= `--min_alt`) or `binom` (alt count significantly above `--seq_error_rate`) |
| `--seq_error_rate` | 0.001 | Per-base error under the null, `--alt_rule binom` only |
| `--alt_binom_alpha` | 0.05 | Significance threshold, `--alt_rule binom` only |
| `--min_alt_floor` | 2 | Never accept a single alt read, `--alt_rule binom` only |
| `--remove_100pct_editing` | true | Remove sites with 100% editing level |
| `--mismatch_window_filter` | false | Remove sites with another mismatch type within w1 nt on either side; w1 is chosen from `--adaptive_windows` to maximise the A-to-G share |
| `--adaptive_windows` | 0,20,40,...,200 | Candidate distances (nt) for the removal of sites near other mismatch types |
| `--clustered_filter` | false | Keep only sites with another site within w2 nt (searched jointly with w1) |
| `--save_intermediates` | 0 | Save the site table of every step from step N on (0 = off; large files) |
| `--requant_editing` | true | Add pileup-based columns (`edit_level_pileup`, `cov_pileup`, `edited_reads`, `ref_reads`) to the published tables, from the BAM that produced each call; original columns are kept |
| `--requant_min_bq` | 30 | Re-quantification: min base quality |
| `--requant_min_mq` | 30 | Re-quantification: min mapping quality |
| `--requant_max_depth` | 1000000 | Re-quantification: `samtools mpileup -d` |

### Per-cell analyses

| Parameter | Default | Description |
|---|---|---|
| `--percell_analyses` | true | Per-cell split + mpileup / AEI (barcode platforms). Set false for caller-comparison runs (the split is memory-heavy) |
| `--percell_calling` | false | Per-cell discovery: run one caller on every cell of one BAM type |
| `--percell_caller` | null | Required with `--percell_calling`: `bcftools`, `reditools3`, `gatk`, `jacusa2` or `scallele` |
| `--percell_bam_type` | null | Required with `--percell_calling`: a sample-, goodCells- or annotated-level BAM type |
| `--percell_aei` | true | Per-cell AEI with RNAEditingIndexer |
| `--percell_chunk_size` | 100 | Cells per RNAEditingIndexer chunk |
| `--rnaeditingindexer_path`, `--aei_genome`, `--aei_regions_bed`, `--aei_snps_bed`, `--aei_refseq_bed`, `--aei_expression_bed` | null | RNAEditingIndexer install and custom inputs |

### Strand modes

| Mode | Callers | Strand of each site |
|---|---|---|
| `annotation` | all | From the gene annotation (`--strand_annots`) |
| `hybrid` | BCFtools, GATK | Read-strand counts, falling back to the annotation when ambiguous |
| `caller` | REDItools3, JACUSA2, scAllele | The caller's own strand call (scAllele: read-cluster strand) |

Modes that do not apply to a caller are skipped. `caller` is not gated on library strandedness: on unstranded libraries it gives what a user of that caller would get.

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

### Large inputs

Caller outputs of 200 MB or more are filtered in chunks of whole chromosomes: all per-site steps run per chunk, and the merged sites then go through the removal of sites near other mismatch types, which needs all sites at once. Outputs are identical to a single run. Environment variables: `DF_CHUNK=0` disables chunking, `DF_CHUNK_MIN_MB` (200) sets the size threshold, `DF_CHUNK_ROWS` (5e6) the minimum rows per chunk and `DF_CHUNK_WORKERS` (1) the number of chunks run at once.

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

Without a count matrix (Smart-seq2), unannotated cells from `called_barcodes` / `cell_id` enter goodCells without the uniform nFeature/nCount/pct_mt QC (`qc_applied = FALSE` in the summary). Supply `count_matrix` to apply those filters.

`recal` and `recalq20` are the coordinate-deduplicated BAM after SplitNCigarReads and BQSR, called at `--qual` and at 20. `recaldefault` is the same BAM called by GATK with its default settings.

Cell-type tier BAMs use meta id `<sample>__<celltype>` so caller/filter publish paths do not collide across cell types of one sample.

## Deduplication strategies

Coordinate and UMI dedup are **one Nextflow task per chromosome** (plus an `_other` / `_all` bucket), then a gather step. When a BAM pools several GEM wells (`@RG SM`), each well is deduplicated separately inside the chromosome task because cell barcodes are unique only within a well.

| Strategy | Flag | Tool | Use case |
|---|---|---|---|
| Coordinate-based | `--dedup coordinate` | GATK MarkDuplicates (`--BARCODE_TAG` from `--dedup_cell_tag`, default CR) | Works for all protocols; Smart-seq2 per-cell BAMs are deduplicated cell by cell |
| UMI-aware | `--dedup umi` | umi_tools dedup (`--per-cell`, cell tag from `--dedup_cell_tag`, `--umi-tag=UB`) | 10x, Drop-seq, and other UMI protocols |
| None | `--dedup none` | Passthrough | Skip deduplication entirely |

Default `--dedup none,coordinate` runs both strategies so RmDup BAM types are produced alongside the unprocessed BAMs.

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
├── genome/                         # STAR index (if built)
├── {sample_id}/
│   ├── 00_validated/               # Validated, MD-tagged, uniquely mapped BAMs
│   ├── 00_qc/                      # Strandedness reports
│   ├── 01_starsolo/                # STARsolo output (if --remap)
│   ├── 01_merged/                  # Lane-merged sample BAM
│   ├── 02_dedup/                   # Deduplicated BAMs + metrics
│   ├── 02_goodcells/               # goodCells barcode list + QC summary
│   ├── 02_barcode_filtered/        # goodCells / annotated filtered BAMs (barcode platforms)
│   ├── 02_bam_types/               # goodCells / annotated / celltype BAMs (Smart-seq2 selective merge)
│   ├── 02_gatk_preprocess/         # SplitNCigar + BQSR BAMs
│   ├── 03_variants/
│   │   ├── bcftools/{bam_type}/    # bcftools mpileup + SNP blacklist
│   │   ├── reditools3/{bam_type}/  # REDItools3 output
│   │   ├── gatk/{bam_type}/        # GATK HaplotypeCaller + SNP blacklist
│   │   ├── jacusa2/{bam_type}/     # JACUSA2 call-1 output
│   │   ├── scallele/{bam_type}/    # SNV TSV + raw VCF + read-cluster table (opt-in)
│   │   └── cellsnp/{bam_type}/     # Cellsnp-lite bulk pileup TSV (opt-in)
│   ├── 03_celltype_bams/           # Per-celltype BAMs
│   ├── 04_filtered/
│   │   └── {caller}/{bam_type}/{strand_mode}/snp_{off,on}/
│   │                               # Filtered sites (pre100 and final), stepwise counts, metrics
│   ├── 04_percell_bams/            # Per-cell BAMs
│   ├── 05_percell_variants/        # Per-cell pileups
│   ├── 05_percell_aei/             # Per-cell Alu Editing Index
│   └── 06_percell_discovery/       # Per-cell discovery (--percell_calling)
├── summary/
│   ├── all_run_metrics.tsv         # One row per configuration
│   ├── all_stepwise_counts.tsv     # Mismatch counts after each filter step
│   ├── rediportal_overlap.tsv      # Known vs novel sites (--rediportal)
│   ├── barcode_coverage_summary.tsv
│   └── barcode_qc_summary.tsv
├── pipeline_trace.txt              # Task-level resource usage
├── pipeline_timeline.html          # Interactive timeline
├── pipeline_report.html            # Execution report
└── pipeline_dag.html               # Workflow DAG
```

## Execution profiles

```bash
-profile singularity # Singularity container (recommended on clusters)
-profile apptainer   # Apptainer container
-profile docker      # Docker container
-profile conda       # Conda environments (see Installation)
-profile standard    # Local execution
-profile sge         # Sun Grid Engine
-profile pbs         # PBS/Torque
-profile slurm       # SLURM
-profile test        # Small resource limits
```

Combine one software profile with one executor profile, e.g. `-profile singularity,slurm`.

Resources are set per process label in `nextflow.config`. All callers share one allocation (12 cpus, 48 GB) so their runtimes are comparable; on SGE, memory is requested with `-l mem` (see the comments in `nextflow.config`).

## Resource tracing

Nextflow's built-in tracing records CPU, memory and wall time for every task in `pipeline_trace.txt` (peak RSS, peak virtual memory, read/write bytes, duration). For callers that run one task per chromosome (scAllele, and JACUSA2 with `--jacusa2_scatter`), the runtime of a configuration is the sum over its chromosome tasks.

## Software versions

Versions used for the benchmark, all in the container and the lock files:

| Software | Version |
|---|---|
| Nextflow | 26.04.6 (>= 23.04 required) |
| samtools | 1.23.1 (htslib 1.23.1) |
| BCFtools | 1.24 |
| GATK | 4.6.2.0 |
| JACUSA2 | 2.1.17 |
| REDItools3 | 3.7 |
| scAllele | 0.0.9.7 |
| Cellsnp-lite | 1.2.3 (htslib 1.24) |
| STAR | 2.7.11b |
| umi_tools | 1.1.6 |
| sinto | 0.9.0 |
| RSeQC | 5.0.5 |
| Python | 3.11.16 (scAllele: 3.10.21) |
| R | 4.5.2; data.table 1.18.6.1, GenomicRanges 1.62.1, Seurat 5.5.1, ggplot2 4.0.3 |

The complete package lists are `envs/*.linux-64.lock.txt`, `envs/scAllele.pip.txt` and `envs/scRNAEditing.r_packages.tsv`.

## Differential editing (separate)

Differential editing analysis (REDIT-LLR, Wilcoxon and limma-based tests with a union consensus) is not part of this pipeline; the scripts and a demonstration are in the code archive of the original study, [doi:10.5281/zenodo.19402220](https://doi.org/10.5281/zenodo.19402220).
