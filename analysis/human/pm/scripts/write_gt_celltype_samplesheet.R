#!/usr/bin/env Rscript
# Build the COLLISION-SAFE per-cell-type samplesheet for the teratoma GT run.
#
# Cell-type discovery must NEVER split a merged multi-teratoma BAM by raw barcode
# (10x barcodes are reused across samples/lanes). Instead we feed the pipeline the
# per-teratoma BAMs (T1-T4) as separate rows under one sample_id = <arm>, each
# with its OWN per-teratoma annotation, and let --split_by_celltype split per
# sample then merge per cell type.
#
# Cell types are filtered on the POOLED cell count across T1-T4 (>= MIN_CELLS);
# per-teratoma annotations are subset to the kept cell types and written to
# gt_barcodes/celltype_flt/. Rows whose filtered annotation is empty are skipped.
#
# NOTE on sample_id here: T1-T4 are DISTINCT biological samples, so grouping them
# under sample_id=<arm> is a DELIBERATE cross-sample POOL for the cell-type tier,
# NOT lanes of one sample. It is only safe because --split_by_celltype splits each
# row (teratoma) by its OWN annotation BEFORE merging per (arm, celltype) — no
# barcode is ever matched across teratomas, so reused 10x barcodes cannot be
# misassigned. Never feed these pooled rows to a barcode-split of a MERGED BAM.
#
# Out: gt_barcodes/celltype_flt/T#_<arm>.tsv  +  samplesheet_pm_gt_celltype.csv

suppressPackageStartupMessages(library(data.table))

base    <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
bcdir   <- file.path(base, "gt_barcodes")
fltdir  <- file.path(bcdir, "celltype_flt")
bamdir  <- file.path(base, "02_gt_bams")
dir.create(fltdir, showWarnings = FALSE, recursive = TRUE)

arms      <- c("ADAR1_KO", "ADARB1_KO", "ADARB2_KO", "AAVS1_ctrl", "NTC_ctrl")
samples   <- c("T1", "T2", "T3", "T4")
MIN_CELLS <- 50

rows <- list()
for (a in arms) {
  # Load all per-teratoma annots for this arm that exist
  perT <- list()
  for (s in samples) {
    f <- file.path(bcdir, sprintf("%s_%s.tsv", s, a))
    if (file.exists(f) && file.size(f) > 0) {
      d <- fread(f, header = FALSE, col.names = c("barcode", "celltype"))
      d[, celltype := gsub(" ", "_", trimws(celltype))]
      perT[[s]] <- d
    }
  }
  if (length(perT) == 0) next

  # Pooled cell count per cell type across teratomas (each barcode = 1 cell/sample)
  pooled <- rbindlist(perT, idcol = "sample")
  cts    <- pooled[, .N, by = celltype][N >= MIN_CELLS, celltype]
  cat(sprintf("== %s: %d cell types >= %d pooled cells ==\n", a, length(cts), MIN_CELLS))
  if (length(cts) == 0) next

  # Write filtered per-teratoma annots (kept cell types only) + samplesheet rows
  for (s in names(perT)) {
    bam <- file.path(bamdir, sprintf("%s_%s.bam", s, a))
    if (!file.exists(bam)) next
    dflt <- perT[[s]][celltype %in% cts]
    if (nrow(dflt) == 0) next
    outf <- file.path(fltdir, sprintf("%s_%s.tsv", s, a))
    fwrite(dflt, outf, sep = "\t", col.names = FALSE)
    rows[[length(rows) + 1]] <- data.table(
      sample_id = a, bam = bam, celltype_annot = outf,
      called_barcodes = outf, platform = "10xv3",
      fastq_r1 = "", fastq_r2 = "")
  }
}

ss <- rbindlist(rows)
fwrite(ss, file.path(base, "samplesheet_pm_gt_celltype.csv"))
cat(sprintf("\nWrote samplesheet_pm_gt_celltype.csv (%d rows across %d arms)\n",
            nrow(ss), length(unique(ss$sample_id))))
print(ss[, .N, by = sample_id])
