# Analysis code

Scripts behind the benchmark figures, the ground-truth definition and the supplementary analyses. The folder layout is the project layout used to produce them, so scripts find each other and their data by relative path:

```
plotting/                     shared ggplot2 theme, colour maps, save helpers; v8 ground-truth paths
publication/<analysis>/scripts/   *_source_data.R / *_plot.R pairs per figure, *_common.R and source-data builders
human/scripts/                ground truth, scoring, ranking and supplementary analyses (R / Python)
human/pm/scripts/             teratoma data preparation and the drivers that run the scripts above
human/<dataset>/scripts/      other datasets (Smart-seq2 / 10x PBMC, DepMap, Replogle)
nextflow -> ..                the pipeline (its bin/ scripts are called by some analyses)
```

`REPO_MANIFEST.tsv` in the project lists every file with its role: `plot` (draws a figure), `source_data` (writes the source data of a figure, or builds a figure-group table), `shared`, `analysis` (upstream analysis or driver) and `data_prep` (download, alignment, guide calling, cell-type annotation, BAM splitting).

## Two scripts per figure

Each figure has two scripts:

- `<name>_source_data.R` computes everything the figure shows from the analysis outputs (tables derived from BAMs, pileups and pipeline calls) and writes it as source data. These scripts are the record of how the source data were made. Most of their inputs are not deposited (raw data, some under restricted access, and BAMs), so they do not run from the Zenodo bundle.
- `<name>_plot.R` reads only that source data and draws the figure. It runs from the Zenodo bundle.

Source data are written to `publication/<analysis>/source_data/<dataset>/<stem>.<part>.tsv`, with one or two tables per panel; large tables are gzipped (`.tsv.gz`). Each table has a `<stem>.<part>.schema.tsv` that records the column types and the order of factor levels, so the plot script gets back exactly the table that was plotted. The helpers `write_sd()` and `read_sd()` are in `plotting/source_data_io.R`. Text fields escape newline, tab and backslash as `\n`, `\t` and `\\`; missing values are `NA`.

New figures follow the same pattern:

1. In `<name>_source_data.R`, source the group's `*_common.R` and `plotting/source_data_io.R`, compute the plotted tables (including any statistics shown on the figure, and any random jitter seed applied to rows in a fixed order), and call `write_sd(table, SRC_DIR, STEM, "<panel>_<what>")` once per table.
2. In `<name>_plot.R`, source the group's `*_plot_common.R`, read each table with `read_part(STEM, part)`, keep the plotting code free of computation beyond layout, and save with `save_figure_sd(plot, STEM, width, height, parts)`. The `.sources.tsv` sidecar then lists the source data tables the figure was drawn from.
3. Add both scripts to `REPO_MANIFEST.tsv` (`source_data` and `plot`) and the source data tables to `ZENODO_FILES.tsv`.

Schematic figures with no data (response figure 14) have a plot script only.

## Regenerating the figures

The source data of every figure are deposited on Zenodo, so the figures can be regenerated without any sequencing data.

```bash
cd analysis
ZENODO_RECORD=<record id> bash download_source_data.sh    # unpacks into this directory and checks MD5 sums
export SCRNA_DIR=$PWD
```

Then run any plot script with the pipeline image (all R packages are in it) or the `scRNAEditing` conda environment (see the main README):

```bash
singularity exec docker://ghcr.io/nytingale/scrnaediting:1.0.0 \
    Rscript publication/teratoma_gt/scripts/TERATOMA_v8_rank_robustness_plot.R
```

Figures are written to `publication/<analysis>/figures/`. The table gives the plot script of each response figure; its source data script has the same name with `_source_data.R` in place of `_plot.R`.

| Response figure | Plot script (`publication/...`) |
|---|---|
| 1 Coverage gates for the ground truth | `gt_definition/scripts/GTDEF_v8_informative_ko_plot.R` |
| 2 Pileup quality, site universe, truth sites per cell type | `gt_definition/scripts/GTDEF_v8_universe_celltype_robustness_plot.R` |
| 3 NTC editing levels of informative negatives | `teratoma_gt/scripts/TERATOMA_v8_negatives_plot.R` |
| 4 Ranks under the nine ground truths | `teratoma_gt/scripts/TERATOMA_v8_rank_robustness_plot.R` |
| 5 Ranks within each truth family | `teratoma_gt/scripts/TERATOMA_v8_rank_robustness_plot.R` |
| 6 Ranking sensitivity to the downstream filters | `teratoma_gt/scripts/TERATOMA_v8_filter_sweep2_rank_rho_plot.R` |
| 7 Rank with and without removal of sites near other mismatch types | `teratoma_gt/scripts/TERATOMA_v8_window_state_rank_scatter_plot.R` |
| 8 Local context filters on candidate sites | `gt_definition/scripts/GTDEF_v8_nomm_cluster_retention_plot.R` |
| 9 Truth-free measures | `teratoma_gt/scripts/TERATOMA_v8_metric_case_plot.R` |
| 10 F1 differences between paired configurations | `teratoma_gt/scripts/TERATOMA_v8_config_factor_plot.R` |
| 11 Editing levels, UMI-aware vs coordinate deduplication | `gt_definition/scripts/GTDEF_hex_dedup_coord_vs_umi_plot.R` |
| 12 Run time, memory and storage per caller | `teratoma_gt/scripts/TERATOMA_v8_caller_resources_plot.R` |
| 13 Ground-truth sites near splice junctions and indels | `gt_definition/scripts/GTDEF_v8_sj_indel_plot.R` |
| 14 Removal of sites near other mismatch types (schematic) | none; taken from Figure 1e of the original manuscript |
| 15 Reference bias correction (WASP) | `wasp_refbias/scripts/WASP_refbias_platforms_plot.R` |
| 16 Full-length vs 3' sites by distance to the 3' end; sites shared between platforms | `three_prime_bias/scripts/DING_3prime_bias_v5_plot.R` |
| 17 JACUSA2 replicate-aware calling | `teratoma_gt/scripts/TERATOMA_v8_replicates_plot.R`, `TERATOMA_v8_replicates_arm_rank_plot.R` |
| 18 Linkage of KO-based sites to SNPs in single cells | `teratoma_gt/scripts/TERATOMA_v8_linkage_snp_plot.R` |
| 19 DepMap CRISPR KO in lung cell lines | `lungTME/scripts/LUNGTME_depmap_exADAR_plot.R` |

The other figures of these analyses follow the same pattern:

| Figures | Plot script (`publication/...`) |
|---|---|
| Ground-truth thresholds (`GTDEF_step1_quality`, `GTDEF_step2_universe`, `GTDEF_robust*`, `GTDEF_panel_sweeps`, `GTDEF_supp_gt_thresholds_v8`, `GTDEF_v8_quality_universe`) | `gt_definition/scripts/GTDEF_gtdefinition_plot.R` |
| Mismatch types and ground-truth sites by BAM processing (`GTDEF_bar_mismatch_spectrum*`, `GTDEF_euler_*`, `GTDEF_upset_*`) | `gt_definition/scripts/GTDEF_bam_processing_plot.R` |
| Benchmark grid (`BENCHMARK_*`), per dataset | `benchmark_grid/scripts/BENCHMARK_grid_plot.R [dataset]` |
| WASP reference bias per platform (`WASP_*`; `jzlxv3` or `jzlss2`) | `wasp_refbias/scripts/WASP_refbias_plot.R [dataset]` |

## Upstream analyses

The `analysis` and `data_prep` scripts start from sequencing data and pipeline outputs, which are not in the Zenodo bundle; they document exactly how every input table was made. The teratoma CRISPR knockout data are GEO GSE248941 (Dailamy et al., 2024). The SGE drivers (`#$` headers, `qsub`) show the order and parameters of each step.

Paths outside the project are read from environment variables whose defaults point inside the pipeline image:

| Variable | Default | Used for |
|---|---|---|
| `SCRNA_DIR` | `.` | project root (this directory) |
| `GT_PILEUP_BIN` | `/opt/conda/envs/gt_pileup/bin` | samtools 1.20 and bedtools 2.31.0, used for the ground-truth pileups |
| `CONDA_ENVS` | `/opt/conda/envs` | the `scRNAEditing`, `scAllele` and `cellsnp` environments |
| `SOFTWARE_DIR` | `/opt/software` | Cell Ranger 10.0.0 and SRA Toolkit 3.0.7 (not in the image; install separately) |

The ground-truth pileups used samtools 1.20 (`gt_pileup` environment), while the pipeline itself uses samtools 1.23.1; both are in the image.
