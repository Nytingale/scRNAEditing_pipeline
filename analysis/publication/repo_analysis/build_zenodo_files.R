#!/usr/bin/env Rscript
# ============================================================================
# build_zenodo_files.R   (2026-11-16)
# Zenodo file list: figure source data and the pipeline inputs we built or curated.
# 1. Figures: each figure's .sources.tsv sidecar (save_figure_sd) plus the .schema.tsv of each table, the tables of the
#    two figures drawn without a sidecar (response 13, 20) and the plot inputs every plot script reads when it loads
#    (cell-type keep lists). FIGS follows the response figures in review/point_by_point_response_2026-11-14.docx; add
#    rows for the main / supplementary figures (paper_figure) when they are fixed. Response figure 18 is a schematic
#    (GTDEF_v8_nomm_schematic_plot.R, no data).
# 2. Pipeline inputs: extended annotations, SNP GRanges, REDIportal VCFs and the STARsolo barcode whitelists
#    (bundle reference_files, unpacked into <bundle_dir>/ by nextflow/references/download_reference_files.sh), and the
#    samplesheets of every pipeline run of the revision with the cell-type and called-barcode tables they name (bundle
#    analysis_data, at their project paths). Not deposited, fetched unchanged by nextflow/references/*/0*_*.sh: genome
#    FASTA, GTF, dbSNP / MGP VCFs, STAR indexes, RSeQC RefSeq BED, main-chromosome lists.
# Output: publication/ZENODO_FILES.tsv
#   path, kind (source_data | schema | plot_input | reference | whitelist | samplesheet | cell_metadata), bytes,
#   scripts, response_figure, paper_figure, env (plot-script environment of the 76-configuration figures),
#   dataset, bundle (analysis_data | reference_files), bundle_dir
#   (read by misc_scripts/build_zenodo_bundle.sh and misc_scripts/build_zenodo_references.sh)
# Usage : Rscript publication/repo_analysis/build_zenodo_files.R   (from anywhere)
# ============================================================================
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("SCRNA_DIR", "."))
ENV76 <- "GT_V8_DIR=human/pm/gt_eval_v9_plot76 TERATOMA_DATASET=pm_gt76 GTDEF_OUT_DATASET=pm_gt76 GT_N_CONFIG=76"
P <- "publication"
FIGS <- fread(text = "
fig|group|dataset|stem|script|env
1|gt_definition|pm_gt|GTDEF_v8_informative_ko|GTDEF_v8_informative_ko_plot|
2|gt_definition|pm_gt76|GTDEF_v8_universe_celltype_robustness|GTDEF_v8_universe_celltype_robustness_plot|76
3|teratoma_gt|pm_gt|TERATOMA_v9_negatives_ntc_editing|TERATOMA_v9_negatives_plot|
4|gt_definition|pm_gt|GTDEF_panel_agtc_v9_response|GTDEF_agtc_v9_plot|
5|gt_definition|pm_gt|GTDEF_bar_truth_spectrum12_response|GTDEF_truth_spectrum12_plot|
6|teratoma_gt|pm_gt76|TERATOMA_v8_rank_robustness_all_tiles|TERATOMA_v8_rank_robustness_plot|76
7|teratoma_gt|pm_gt76|TERATOMA_v8_rank_robustness_by_family|TERATOMA_v8_rank_robustness_plot|76
8|teratoma_gt|pm_gt76|TERATOMA_v8_rank_gap_top5|TERATOMA_v8_rank_gap_plot|76
9|teratoma_gt|pm_gt76|TERATOMA_v8_caller_margin|TERATOMA_v8_caller_margin_plot|76
10|teratoma_gt|pm_gt|TERATOMA_v8_filter_sweep2_rank_rho_pwcov|TERATOMA_v8_filter_sweep2_rank_rho_plot|
11|teratoma_gt|pm_gt76|TERATOMA_v8_window_state_rank_scatter|TERATOMA_v8_window_state_rank_scatter_plot|76
12|gt_definition|pm_gt|GTDEF_v8_nomm_cluster_retention|GTDEF_v8_nomm_cluster_retention_plot|
14|teratoma_gt|pm_gt|TERATOMA_v8_config_factor_delta_f1|TERATOMA_v8_config_factor_plot|
15|gt_definition|pm_gt|GTDEF_hex_editlevel_coord_vs_umi|GTDEF_hex_dedup_coord_vs_umi_plot|
16|teratoma_gt|pm_gt|TERATOMA_v8_caller_resources|TERATOMA_v8_caller_resources_plot|
17|gt_definition|pm_gt|GTDEF_v8_sj_indel|GTDEF_v8_sj_indel_plot|
19|wasp_refbias|combined|WASP_panel_refbias_platforms|WASP_refbias_platforms_plot|
21|teratoma_gt|pm_gt76|TERATOMA_v8_replicates_arm_rank|TERATOMA_v8_replicates_arm_rank_plot|76
22|teratoma_gt|pm_gt|TERATOMA_v8_linkage_snp_cells|TERATOMA_v8_linkage_snp_plot|
23|lungTME|dl|LUNGTME_panel_depmap_adar_isg_exADAR|LUNGTME_depmap_exADAR_plot|
", sep = "|", colClasses = "character")
FIGS[, env := fifelse(env == "76", ENV76, "")]

sd <- rbindlist(lapply(seq_len(nrow(FIGS)), function(i) {
  f <- FIGS[i]
  side <- file.path(P, f$group, "source_data", f$dataset, paste0(f$stem, ".sources.tsv"))
  if (!file.exists(side)) stop("no sidecar: ", side)
  s <- fread(side)
  if (!all(s$source_exists)) stop("sidecar lists missing tables: ", side)
  data.table(path = s$source_file, kind = "source_data", scripts = f$script, response_figure = f$fig, env = f$env)
}))
extra <- rbind(
  data.table(path = file.path(P, "teratoma_gt/source_data/pm_gt76/TRUTHFREE_v10_harmonic_scatter.tsv"),
             kind = "source_data", scripts = "TRUTHFREE_v10_harmonic_scatter", response_figure = "13", env = ""),
  data.table(path = file.path(P, "three_prime_bias/source_data",
                              paste0("ding_3prime_bias_v5.", c("edited_by_region", "dist3p_bins", "overlap_by_region"), ".tsv")),
             kind = "source_data", scripts = "DING_3prime_bias_v5_plot", response_figure = "20", env = ""))
sd <- rbind(sd, extra)
schema_of <- function(p) sub("\\.tsv(\\.gz)?$", ".schema.tsv", p)
sch <- sd[file.exists(schema_of(path))][, `:=`(path = schema_of(path), kind = "schema")]
pin <- data.table(path = c("human/pm/gt_eval_v8/unit_cov/celltype_keep_pwcov.txt",
                           "human/pm/gt_eval_v9_plot76/unit_cov/celltype_keep_pwcov.txt",
                           "human/pm/gt_eval_v9/unit_cov/celltype_keep_pwcov.txt"),
                  kind = "plot_input", scripts = "every TERATOMA / GTDEF plot script (plotting/gt_v8_paths.R V8_CT_KEEP)",
                  response_figure = "", env = c("", ENV76, "target of the gt_eval_v9_plot76 link"))
z <- rbind(sd, sch, pin)
z <- z[, .(kind = kind[1], scripts = paste(unique(unlist(strsplit(scripts, ","))), collapse = ","),
           response_figure = paste(unique(response_figure[nzchar(response_figure)]), collapse = ","),
           env = paste(unique(env[nzchar(env)]), collapse = "; ")), by = path]
z[, `:=`(paper_figure = "", dataset = "", bundle = "analysis_data", bundle_dir = "")]
z[, o := suppressWarnings(as.integer(sub(",.*", "", response_figure)))]
setorder(z, o, kind, path, na.last = TRUE); z[, o := NULL]

# ── Pipeline inputs ──────────────────────────────────────────────────────────
GENOMES <- c("GRCh38_analysis_set", "mm39", "dm6", "danRer11", "xenTro10", "xenLae10", "xcOctVulg1.2", "obim")
ref <- rbind(
  data.table(path = sprintf("annotations/annotation_files/extended_annots_%s.RData", GENOMES), kind = "reference",
             scripts = "--strand_annots / --alu_annots", bundle_dir = "annotation_files"),
  data.table(path = c("annotations/SNP_files/dbSNP156_common_GRanges.RData", "annotations/SNP_files/GRCm39_snps_GRanges.RData",
                      sprintf("annotations/SNP_files/%s_SNPs_GRanges.RData", GENOMES[-(1:2)])),
             kind = "reference", scripts = "--snp_file", bundle_dir = "SNP_files"),
  data.table(path = c("human/ref_annots/REDIportal_hg38_v2.vcf.gz", "human/ref_annots/REDIportal_hg38_v2.vcf.gz.tbi",
                      "mouse/ref_annots/REDIportal_mm39.vcf"), kind = "reference", scripts = "--rediportal", bundle_dir = "rediportal"),
  data.table(path = file.path("nextflow/assets/whitelists",
                              c("737K-august-2016.txt", "3M-february-2018.txt", "inDropCB1.txt", "inDropCB2.txt", "README.txt")),
             kind = "whitelist", scripts = "STARsolo --soloCBwhitelist (--remap true)", bundle_dir = "whitelists"))
ref[, `:=`(dataset = fcase(grepl("GRCh38|dbSNP|hg38", path), "GRCh38_analysis_set", grepl("mm39|GRCm39", path), "mm39",
                           kind == "whitelist", "",
                           default = sub("^extended_annots_|_SNPs_GRanges$", "", sub("\\.RData$", "", basename(path)))),
           bundle = "reference_files")]

# Samplesheets of the revision runs (run_pipeline_pm_gt_v8.sh, run_pipeline_aorta_v8.sh, run_pipeline_jzl*.sh) and the
# cell-type / called-barcode tables they name. BAM paths in them are project paths; the BAMs are not deposited.
SHEETS <- c(
  sprintf("human/pm/samplesheet_pm_gt_%s.csv", c("merged_dedupv6", "persample_dedupv6", "celltype_as_sample_dedupv6",
                                                 "merged_dedupv6_umimd_ntc", "persample_dedupv6_umimd_ntc",
                                                 "celltype_as_sample_dedupv6_umimd_ntc", "persample", "celltype_v8_NTC")),
  sprintf("human/pm/samplesheet_umimd_%s.csv", c("arm", "arm_ctrl", "sample", "celltype")),
  sprintf("human/%s/samplesheet_%s.csv", c("jzlxv3", "jzlss2", "jzl10x", "jzlcel"), c("jzlxv3", "jzlss2", "jzl10x", "jzlcel")),
  sprintf("mouse/aorta_gse255661/pipeline_bams/samplesheet_%s.csv", c("aorta_ctrl", "aorta_ctrl_umimd")),
  sprintf("mouse/aorta_gse255661/samplesheet_%s.csv", c("umimd_aorta", "umimd_aorta_ctrl")))
if (!all(file.exists(SHEETS))) stop("missing samplesheets: ", paste(SHEETS[!file.exists(SHEETS)], collapse = ", "))
# The deposited samplesheets are copies with project-relative paths, in publication/samplesheets/<dataset>/.
ROOT <- paste0(getwd(), "/")
SS_DIR <- file.path(P, "samplesheets")
cm <- rbindlist(lapply(SHEETS, function(s) {
  x <- fread(s, colClasses = "character")
  ds <- sub("^(human|mouse)/([^/]+)/.*", "\\2", s)
  for (k in names(x)) set(x, j = k, value = gsub(ROOT, "", x[[k]], fixed = TRUE))
  abs <- unlist(lapply(x, function(v) grep("^/", v, value = TRUE)))
  if (length(abs)) stop(s, ": absolute paths outside the project: ", paste(head(abs, 3), collapse = ", "))
  out <- file.path(SS_DIR, ds, basename(s))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, out)
  cols <- intersect(names(x), c("celltype_annot", "called_barcodes", "whitelist_files"))
  f <- unique(unlist(lapply(cols, function(k) unlist(strsplit(x[[k]], ";")))))
  f <- f[nzchar(f)]
  rbind(data.table(path = out, kind = "samplesheet"), data.table(path = f, kind = rep("cell_metadata", length(f))))[
    , dataset := ds]
}))
cm <- cm[, .(kind = kind[1], dataset = paste(unique(dataset), collapse = ",")), by = path]
cm[, `:=`(scripts = fifelse(kind == "samplesheet", "--input", "samplesheet celltype_annot / called_barcodes"),
          bundle = "analysis_data", bundle_dir = "")]
if (cm[kind == "cell_metadata", any(startsWith(path, "/"))]) stop("cell metadata outside the project: ",
                                                                   paste(cm[startsWith(path, "/"), path], collapse = ", "))
z <- rbind(z, ref, cm, fill = TRUE)
for (k in names(z)) if (is.character(z[[k]])) set(z, which(is.na(z[[k]])), k, "")

miss <- z[!file.exists(path), path]
if (length(miss)) stop("missing files:\n", paste(miss, collapse = "\n"))
if (anyDuplicated(z$path)) stop("duplicated paths: ", paste(z$path[duplicated(z$path)], collapse = ", "))
z[, bytes := file.size(path)]
setcolorder(z, c("path", "kind", "bytes", "scripts", "response_figure", "paper_figure", "env", "dataset", "bundle", "bundle_dir"))
fwrite(z, file.path(P, "ZENODO_FILES.tsv"), sep = "\t", quote = FALSE)
cat(sprintf("%d response figures; %d files, %.1f MB in total\n",
            uniqueN(unlist(strsplit(z[nzchar(response_figure), response_figure], ","))), nrow(z), sum(z$bytes) / 1024^2))
print(z[, .(files = .N, MB = round(sum(bytes) / 1024^2, 1)), by = .(bundle, kind)])
