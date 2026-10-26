#!/usr/bin/env Rscript
# Write per-sample × arm barcode TSVs and Nextflow samplesheets for the
# paper-faithful teratoma ground-truth contrast.
#
# Inputs:  seurat_analysis/ground_truth_cells_paper.tsv
#          01_starsolobams/T{1-4}_merged.bam  (split later; paths recorded now)
# Outputs: gt_barcodes/{sample}_{arm}.tsv     barcode \t celltype
#          samplesheet_pm_gt_persample.csv    20 rows (T1–T4 × 5 arms)
#          samplesheet_pm_gt_merged.csv       5 rows  (arms merged across T1–T4)
#
# Merged BAMs are for pseudobulk depth (ADAR1_KO 680 cells / 5.2M UMI).
# Per-sample BAMs are kept as JACUSA2/limma replicates. Controls stay separate
# (AAVS1_ctrl vs NTC_ctrl). Do not use the old monolithic T1–T4 samplesheet.
#
# LANE vs SAMPLE (pipeline contract): in a samplesheet, sample_id is the TRUE
# biological sample and the sample-level analysis partition. Multiple rows sharing
# a sample_id are LANES of that SAME sample (same cells -> same barcodes), merged
# before analysis. DISTINCT biological samples MUST get DISTINCT sample_ids: 10x
# barcodes are reused across samples, so merging different samples under one id and
# doing a barcode-based split would MISASSIGN reads. That is exactly why:
#   - persample sheet = 20 true partitions (T1-T4 x 5 arms), one sample_id each; and
#   - cross-sample cell-type pooling (T1-T4 -> arm) is NOT done here by faking a
#     merged sample_id, but via the celltype sheet + pipeline --split_by_celltype
#     (split each teratoma by cell type BEFORE merge; see write_gt_celltype_samplesheet.R).

suppressPackageStartupMessages(library(data.table))

base   <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
gt     <- fread(file.path(base, "seurat_analysis/ground_truth_cells_paper.tsv"))
bamdir <- file.path(base, "01_starsolobams")
splitdir <- file.path(base, "02_gt_bams")
bcdir  <- file.path(base, "gt_barcodes")
dir.create(bcdir, showWarnings = FALSE, recursive = TRUE)
dir.create(splitdir, showWarnings = FALSE, recursive = TRUE)

arms <- c("ADAR1_KO", "ADARB1_KO", "ADARB2_KO", "AAVS1_ctrl", "NTC_ctrl")
samples <- c("T1", "T2", "T3", "T4")
gt <- gt[arm %in% arms]
stopifnot(nrow(gt) > 0)

# barcode \t celltype  (SINTO_EXTRACT_QC uses col1; celltype split uses both)
for (s in samples) {
  for (a in arms) {
    d <- gt[sample == s & arm == a, .(barcode, celltype)]
    fwrite(d, file.path(bcdir, sprintf("%s_%s.tsv", s, a)), sep = "\t", col.names = FALSE)
  }
}

# Merged annot: prefix barcode with sample so CBs are unique across GEM wells
# (needed if the merged BAM is later barcode-filtered). Split BAMs keep raw CBs.
merged_annot <- gt[, .(barcode = paste(sample, barcode, sep = "_"), celltype)]
for (a in arms) {
  fwrite(gt[arm == a, .(barcode, celltype)],
         file.path(bcdir, sprintf("merged_%s.tsv", a)), sep = "\t", col.names = FALSE)
}

mkrow <- function(id, bam, annot, called) {
  data.table(sample_id = id, bam = bam, celltype_annot = annot,
             called_barcodes = called, platform = "10xv3",
             fastq_r1 = "", fastq_r2 = "")
}

persample <- rbindlist(lapply(samples, function(s) {
  rbindlist(lapply(arms, function(a) {
    bc <- file.path(bcdir, sprintf("%s_%s.tsv", s, a))
    bam <- file.path(splitdir, sprintf("%s_%s.bam", s, a))
    mkrow(paste(s, a, sep = "_"), bam, bc, bc)
  }))
}))
fwrite(persample, file.path(base, "samplesheet_pm_gt_persample.csv"))

merged <- rbindlist(lapply(arms, function(a) {
  bc <- file.path(bcdir, sprintf("merged_%s.tsv", a))
  bam <- file.path(splitdir, sprintf("merged_%s.bam", a))
  mkrow(a, bam, bc, bc)
}))
fwrite(merged, file.path(base, "samplesheet_pm_gt_merged.csv"))

cat("=== cells per sample x arm ===\n")
print(dcast(gt[, .N, by = .(arm, sample)], arm ~ sample, value.var = "N", fill = 0))
cat("\nWrote", length(list.files(bcdir, "*.tsv")), "barcode TSVs in", bcdir, "\n")
cat("Wrote samplesheet_pm_gt_persample.csv (", nrow(persample), "rows )\n")
cat("Wrote samplesheet_pm_gt_merged.csv    (", nrow(merged), "rows )\n")
