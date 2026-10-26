# Reference files

How every reference and annotation file given to the pipeline was obtained or built. Files that are downloaded unchanged (genome, GTF, dbSNP VCF, REDIportal table, RSeQC BED) are fetched by the scripts here. Files we built (`extended_annots_*.RData`, `*_SNPs_GRanges.RData`) can be rebuilt with these scripts, or downloaded ready-made from Zenodo (see [Download the built files](#download-the-built-files)).

## Overview (human, GRCh38 analysis set)

| Parameter | File | Source | Built by |
|---|---|---|---|
| `--genome_fasta` | `hg38.analysisSet.fa` (+ `.fai`, `.dict`) | [UCSC hg38 analysis set](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/) (no ALT contigs) | `human/01_download_genome_and_gtf.sh` |
| `--genome_gtf` | `hg38.ncbiRefSeq.gtf` | [UCSC NCBI RefSeq track](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/) | `human/01_download_genome_and_gtf.sh` |
| `--main_chroms` | `chrlist_main_chrs.txt` (chr1-22, X, Y) | from the FASTA index | `human/01_download_genome_and_gtf.sh` |
| `--known_sites_vcf` | `dbSNP156_GRCh38.vcf.gz` | [dbSNP build 156](https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/), RefSeq accessions renamed to UCSC names with the NCBI assembly report | `human/02_download_dbsnp_vcf.sh` |
| `--star_index` | `stardb_hg38_<overhang>/` | STAR `genomeGenerate` from the FASTA and GTF | `human/03_build_star_index.sh` |
| `--rediportal` | `REDIportal_hg38_v2.vcf.gz` | [REDIportal v2](http://srv00.recas.ba.infn.it/atlas/) TABLE1 (hg38) | `human/04_rediportal_vcf.sh` |
| `--strand_annots`, `--alu_annots` | `extended_annots_GRCh38_analysis_set.RData` | [NCBI RefSeq GTF of the GRCh38 analysis set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) + UCSC RepeatMasker (`rmsk.txt.gz`) | `human/05_alu_from_rmsk.sh`, `human/06_extended_annots.sh` |
| `--refseq_bed` | `hg38_RefSeq.bed` | [RSeQC](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/) BED12 | `human/07_refseq_bed.sh` |
| `--snp_file` | `dbSNP156_common_GRanges.RData` | common SNPs of dbSNP build 156 | `human/08_dbsnp_common.sh` |

Mouse (mm39) uses the same parameters with `mouse/` scripts (see [Mouse](#mouse-mm39)).

## Software

Downloads, renaming, indexing and the R steps use the pipeline's main environment (`envs/scRNAEditing.*`, or the container): wget, samtools, bcftools, tabix / bgzip, GATK, STAR, R with GenomicRanges, rtracklayer and data.table. Only the extended annotations also need `gtf2bed` (BEDOPS) and bedtools, in the small `refbuild` environment:

```bash
mamba create -n refbuild --file envs/refbuild.linux-64.lock.txt    # BEDOPS 2.4.42, bedtools 2.31.0
# or, in the container:  singularity exec <image> micromamba run -p /opt/conda/envs/refbuild ...
```

Every script writes into `$REFDIR` (default: the current directory) and skips downloads that are already there.

## Human (GRCh38 analysis set)

```bash
export REFDIR=$PWD/ref_hg38
bash references/human/01_download_genome_and_gtf.sh   # FASTA (3.1 GB), .fai, .dict, GTF, chromosome list
bash references/human/02_download_dbsnp_vcf.sh        # dbSNP 156 VCF (~25 GB), chr names -> UCSC; several hours
bash references/human/03_build_star_index.sh 100      # ~32 GB RAM, 12 threads; also 49 for 50-nt reads
bash references/human/04_rediportal_vcf.sh            # REDIportal v2 -> VCF (15,681,871 sites)
bash references/human/05_alu_from_rmsk.sh             # Alu elements per strand from RepeatMasker
bash references/human/06_extended_annots.sh           # needs refbuild on PATH -> extended_annots_GRCh38_analysis_set.RData
bash references/human/07_refseq_bed.sh                # RSeQC hg38_RefSeq.bed
bash references/human/08_dbsnp_common.sh              # dbSNP156_common_GRanges.RData
```

**Extended gene annotations** (`--strand_annots`). The pipeline strands each site by the gene annotation and, where annotations on both strands overlap it, by the region type on each strand (exonic features, i.e. exon, UTR, CDS, start / stop codon, ranked equally, over intron over intergenic; ties stay unstranded). The human annotation is built from NCBI's RefSeq GTF for the GRCh38 analysis set (`GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz`, UCSC chromosome names; it has gene rows and no UTR rows), not from the UCSC `hg38.ncbiRefSeq.gtf` used for alignment; mouse uses the UCSC mm39 `ncbiRefSeq` GTF. This needs every position of the genome assigned on each strand, so `common/generate_extended_bed.sh` converts the GTF to BED (`gtf2bed`), keeps the main chromosomes (drops `_fix`, `_alt`, `_random`, `chrUn`, `chrM`), splits it by strand and adds, per strand, introns (transcript bodies minus merged exons, `bedtools subtract`) and intergenic regions (the genome not covered by a transcript, `bedtools complement`). `common/make_extended_annots.R` reads the two BEDs into GRanges `plus.annot` and `minus.annot` (1-based, columns `id` = gene name and `geneType` = the GTF feature type, i.e. exon, CDS, start / stop codon, transcript, plus gene (NCBI GTF) or 5UTR / 3UTR (UCSC GTF), or the added intron or intergenic; the strand is stored as `*` because the object name gives it).

**Alu annotations** (`--alu_annots`, human only). RepeatMasker rows containing "Alu", split by strand, as `alu.plus.annot` and `alu.minus.annot` (column `id` = repeat name) in the same RData. Sites still unstranded after the gene annotation get the strand of the Alu element they lie in, if Alu elements overlap them on one strand only.

**Known SNPs** (`--snp_file`). An RData file with a GRanges `snp_gr` of width-1 positions; sites at these positions are removed (`--snp_subtract_modes on`). Human: the common variants of dbSNP build 156, i.e. every record of the dbSNP 156 VCF with `INFO/COMMON=1`, all variant classes (20,290,286 SNVs, 3,292,079 insertions and deletions and 156 MNVs; 23,582,521 records at 22,993,765 unique positions; `08_dbsnp_common.sh` reproduces the original study's `dbSNP156_common.RData` exactly, in the same order). Mouse: every SNP position of the Mouse Genomes Project REL-2112-v8 VCF on the main chromosomes (83,212,162; identical to the original study's `mm39_SNPs.RData`).
`common/snp_list_to_granges.R` converts the original study's SNP lists (data frame `snplist`, column `site` = `chr:pos`) to this format; the pipeline also accepts those lists and a BED file directly.

**GATK recalibration** (`--known_sites_vcf`) uses the full dbSNP 156 VCF, not the common subset.

## Mouse (mm39)

```bash
export REFDIR=$PWD/ref_mm39
bash references/mouse/01_download_genome_and_gtf.sh   # UCSC mm39 FASTA, ncbiRefSeq GTF, chrom.sizes, chr1-19, X, Y
bash references/mouse/02_download_mgp_snp_vcf.sh      # Mouse Genomes Project REL-2112-v8 SNPs (--known_sites_vcf), chr names -> UCSC
bash references/mouse/03_make_snp_rdata.sh            # GRCm39_snps_GRanges.RData (--snp_file): all MGP SNP positions on the main chromosomes
bash references/mouse/04_build_star_index.sh 49
bash references/mouse/05_extended_annots.sh           # needs refbuild on PATH -> extended_annots_mm39.RData (no Alu)
bash references/mouse/06_refseq_bed.sh                # mm39_RefSeq.bed (BED12 from the UCSC ncbiRefSeq table)
Rscript references/mouse/07_rediportal_mm39.R         # REDIportal mouse (mm10, ExperimentHub EH8233) lifted to mm39
```

## Other genomes

The original study ([doi:10.5281/zenodo.19402220](https://doi.org/10.5281/zenodo.19402220)) also provided extended annotations and SNP lists for dm6, danRer11, xenTro10, xenLae10, xcOctVulg1.2 (*Octopus vulgaris*) and ASM119413v2 (*Octopus bimaculoides*, named `obim`). The extended annotations were built with the same two scripts; only the annotation source and the main chromosomes differ. The SNP lists were built in 2024 as described in the original study's Methods; those steps are not scripted here, and the built files are on Zenodo.

| Genome | Genome and gene annotation | Main chromosomes | SNP list |
|---|---|---|---|
| dm6 | UCSC (ncbiRefSeq) | chr2L, chr2R, chr3L, chr3R, chr4, chrX, chrY | *Drosophila melanogaster* Genetic Reference Panel 2 (accessed 2024-05-21), lifted from dm3 to dm6 with UCSC liftOver |
| danRer11 | UCSC (ncbiRefSeq) | chr1-25 | variants of all strains in SNPfisher, lifted from danRer7 to danRer11 (UCSC chain), 1:1 single-nucleotide conversions only |
| xenTro10 | Xenbase | chr1-10 | called from whole-genome sequencing (SRR5110206-SRR5110208) |
| xenLae10 | Xenbase | chr1-8 L / S, chr9_10L, chr9_10S | called from whole-genome sequencing (PRJNA313213) |
| xcOctVulg1.2 | NCBI (submitter annotation; files named `embl`) | chr1-30 | called from whole-genome sequencing (ERR11272671-ERR11272678) |
| obim (ASM119413v2) | NCBI RefSeq | chr1-30 | called from whole-genome sequencing (PRJNA270931, PRJEB51312, PRJNA790947) |

Whole-genome calls: reads trimmed with fastp 0.23.4, aligned with STAR 2.7.11a, duplicates marked with GATK MarkDuplicates, variants called with FreeBayes 1.3.6; single-nucleotide variants with quality ≥ 20 and depth ≥ 10 kept.

For a new genome: put `<prefix>.<annot>.gtf.gz` and `<prefix>.chrom.sizes` in a directory, run `common/generate_extended_bed.sh <prefix> <annot>` there, then `common/make_extended_annots.R` with a main-chromosome list.

## Download the built files

`bash references/download_reference_files.sh <dir>` downloads `scRNAEditing_reference_files_1.0.0.tar.gz` (452 MB) from Zenodo (record: ZENODO_RECORD_TBD), checks its MD5 sums and unpacks `annotation_files/extended_annots_<genome>.RData` and `SNP_files/<genome>_SNPs_GRanges.RData` (human `dbSNP156_common_GRanges.RData`, mouse `GRCm39_snps_GRanges.RData`) for all eight genomes, and `rediportal/REDIportal_hg38_v2.vcf.gz` (+ `.tbi`) and `rediportal/REDIportal_mm39.vcf` for `--rediportal`. These are the exact files used for the paper. The other public inputs (genome, GTF, dbSNP VCF, RSeQC BED) come from the download scripts above.

## Checks

Rebuilt on 2026-10-24 with these scripts and compared with the files used for the paper:

| File | Result |
|---|---|
| `REDIportal_hg38_v2.vcf` | byte-identical |
| `hg38_RefSeq.bed` | byte-identical (RSeQC file, md5 `8c83e959…`) |
| `mm39_RefSeq.bed` | byte-identical |
| `chrlist_main_chrs.txt` (human) | byte-identical |
| `rediportal_full_mm10.rda` (ExperimentHub download) | byte-identical |
| `extended_annots_mm39.RData` | same ranges and annotations, in the same order |
| `extended_annots_GRCh38_analysis_set.RData`, Alu objects | same ranges and annotations, in the same order |
| `extended_annots_GRCh38_analysis_set.RData`, gene objects | not identical: NCBI has replaced the GTF since 2024 (0.2% of the 2024 rows absent, 0.5% added, from annotation updates); use the Zenodo file to reproduce the paper |
| `dbSNP156_common.RData` and `dbSNP156_common_GRanges.RData` from the dbSNP 156 VCF (`08_dbsnp_common.sh`) | identical, same order; GRanges file byte-identical |
| `*_SNPs_GRanges.RData` (dbSNP156_common, dm6, danRer11) from the original lists | identical |

`convert2bed` (BEDOPS 2.4.39 to 2.4.42) can abort with heap corruption after writing all records; `generate_extended_bed.sh` therefore checks the record count instead of the exit status.
