#!/usr/bin/env Rscript
# teratoma_annotation.R — Seurat preprocessing + label transfer for ADAR-KO teratomas
# Dailamy et al. 2024; WT reference McDonald, Wu et al. 2020 (GSE156170)
#
# Biological steps follow the published GitHub scripts (Seurat 5 API names only):
#   Phase D/F: SammiLyu/scScreens_ADARs Figure6/F06_01_ADARKO_Preprocessing.R
#   Phase E:   yanwu2014/teratoma-analysis-code Figure1/01_human_clustering.R
#              plus published cluster labels from Figure1/dm-ter-human-merged_metadata.tsv
#              (F06_01 maps onto WT4H1_YWsubset_integratedObj; we rebuild that
#              integration then attach those labels).
#
# Seurat 5 / cluster runtime only (not biology):
#   - Create then subset so CreateSeuratObject cannot drop all cells (LogMap)
#   - Harmonize WT barcodes "-" vs "." so published clusters match
#   - Checkpoint D/E/F; GEM species call (UMI fallback only if file is absent)
#   - ncol()/nrow() are double in Seurat 5; n_cells()/n_features() are integer
#   - future 1.75: do not cap future.globals.maxSize (MapQuery ~16.1 GiB)
#   - DoubletFinder: current paramSweep/doubletFinder names replace the
#     published _v3 names; retain the published sct=TRUE and PCs=1:75.
#
# Phases:
#   D — QC, DoubletFinder, human singlets
#   E — WT reference (integration + published cluster labels)
#   F — Label transfer
#   G — barcode → celltype TSVs for the pipeline
# Plots live in plotting/plot_teratoma_annotation.R (this script writes TSVs only).
#
# Resume: skip a finished phase unless FORCE=1
# Threads: NSLOTS or ANNOT_THREADS (default 4)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DoubletFinder)
  library(data.table)
  library(Matrix)
})

n_threads <- as.integer(Sys.getenv("NSLOTS", Sys.getenv("ANNOT_THREADS", "4")))
if (is.na(n_threads) || n_threads < 1) n_threads <- 1
force_rerun <- identical(Sys.getenv("FORCE", ""), "1")

# future 1.75: setting future.globals.maxSize puts a cap on sequential
# plans too (maxSizeOfObjects). MapQuery exports ~16.1 GiB. Leave
# sequential uncapped; do not set a 16 GiB option.
configure_seurat_future <- function() {
  options(future.globals.maxSize = Inf)
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential, maxSizeOfObjects = Inf)
    log_msg("future plan=sequential  maxSizeOfObjects=Inf")
  }
}

cr_dir  <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/cellranger_out")
wt_dir  <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/wt_reference")
out_dir <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/seurat_analysis")
pm_dir  <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples  <- c("T1", "T2", "T3", "T4")
ter_list <- c("ter1", "ter2", "ter3", "ter4")

human_prefix <- "GRCh38-"
mouse_prefix <- "GRCm38-"

phase_d_rda <- file.path(out_dir, "KO_4ter_DFed_humanCells_list.rda")
phase_e_rds <- file.path(out_dir, "wt_reference.rds")
phase_f_rda <- file.path(out_dir, "KO_4ter_annotated_list.rda")

log_msg <- function(...) cat(paste0(...), "\n")
configure_seurat_future()

# Seurat 5 / SeuratObject 5 return double from ncol()/nrow().
n_cells <- function(obj) as.integer(ncol(obj))
n_features <- function(obj) as.integer(nrow(obj))

require_file <- function(path, what) {
  if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)
  path
}

require_ncells <- function(obj, label, min_n = 1L) {
  n <- n_cells(obj)
  if (is.null(n) || is.na(n) || n < min_n) {
    stop(label, " has ", if (is.null(n)) "NULL" else n, " cells (need >= ", min_n, ")",
         call. = FALSE)
  }
  invisible(n)
}

sanitize_celltype <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Unknown"
  x <- gsub("[/ ]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

as_named_sparse <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    stop("count matrix has missing rownames", call. = FALSE)
  }
  rownames(mat) <- make.unique(as.character(rownames(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  colnames(mat) <- make.unique(as.character(colnames(mat)))
  mat
}

# Seurat 5.5 LogMap dies if CreateSeuratObject's min.features/min.cells
# drops every cell or every gene. Always create first, subset after.
# F06_01 / yanwu CreateSeuratObject(min.cells = round(0.001 * n))
filter_min_cells <- function(obj, frac = 0.001) {
  min_c <- max(1L, round(frac * n_cells(obj)))
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  keep <- Matrix::rowSums(counts > 0) >= min_c
  n_keep <- sum(keep)
  if (n_keep < 100L) {
    stop("min.cells=", min_c, " left ", n_keep, " genes", call. = FALSE)
  }
  subset(obj, features = rownames(obj)[keep])
}

new_seurat <- function(counts, meta = NULL, project = "SeuratProject") {
  counts <- as_named_sparse(counts)
  if (!is.null(meta)) {
    meta <- as.data.frame(meta)
    miss <- setdiff(colnames(counts), rownames(meta))
    if (length(miss) > 0) {
      extra <- data.frame(row.names = miss)
      meta <- rbind(meta, extra)
    }
    meta <- meta[colnames(counts), , drop = FALSE]
  }
  CreateSeuratObject(
    counts = counts, meta.data = meta, project = project,
    min.cells = 0, min.features = 0
  )
}

extract_counts <- function(raw) {
  if (is.list(raw) && !is.data.frame(raw)) {
    if ("Gene Expression" %in% names(raw)) return(raw[["Gene Expression"]])
    return(raw[[1]])
  }
  raw
}

n_pcs_available <- function(obj) {
  if (!"pca" %in% Reductions(obj)) return(0L)
  as.integer(ncol(Embeddings(obj, reduction = "pca")))
}

parse_cr_metric <- function(metrics_csv, colname) {
  m <- fread(metrics_csv)
  if (!colname %in% names(m)) return(NA_real_)
  val <- as.character(m[[colname]][1])
  val <- gsub(",", "", val)
  val <- gsub("%", "", val)
  suppressWarnings(as.numeric(val))
}

preflight_cellranger <- function() {
  log_msg("Preflight: CellRanger outs")
  for (s in samples) {
    mtx_dir <- file.path(cr_dir, s, "outs", "filtered_feature_bc_matrix")
    mtx <- file.path(mtx_dir, "matrix.mtx.gz")
    if (!file.exists(mtx)) mtx <- file.path(mtx_dir, "matrix.mtx")
    require_file(mtx, paste(s, "filtered matrix"))
    metrics <- file.path(cr_dir, s, "outs", "metrics_summary.csv")
    require_file(metrics, paste(s, "metrics_summary.csv"))
    n_cells <- parse_cr_metric(metrics, "Estimated Number of Cells")
    med_genes <- parse_cr_metric(metrics, "Median Genes per Cell")
    log_msg(sprintf("  %s: estimated_cells=%s  median_genes=%s",
                    s, n_cells, med_genes))
    if (!is.na(med_genes) && med_genes < 50) {
      stop(s, " median genes/cell is ", med_genes,
           " — CellRanger looks empty (wrong pairing or failed count)",
           call. = FALSE)
    }
    if (!is.na(n_cells) && n_cells < 200) {
      stop(s, " estimated cells is ", n_cells,
           " — too few for DoubletFinder / label transfer",
           call. = FALSE)
    }
  }
}

read_gem_calls <- function(sample, ter_name, cell_names) {
  candidates <- c(
    file.path(cr_dir, sample, "outs", "analysis", "gem_classification.csv"),
    file.path(cr_dir, sample, "outs", "gem_classification.csv")
  )
  gem_csv <- candidates[file.exists(candidates)][1]
  if (is.na(gem_csv)) return(NULL)
  gem <- fread(gem_csv)
  bc_col <- intersect(c("barcode", "Barcode"), names(gem))[1]
  call_col <- intersect(c("call", "Call"), names(gem))[1]
  if (is.na(bc_col) || is.na(call_col)) return(NULL)
  keys <- paste0(ter_name, "_", gem[[bc_col]])
  calls <- as.character(gem[[call_col]])
  names(calls) <- keys
  calls[intersect(cell_names, names(calls))]
}

human_frac_calls <- function(obj) {
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  h_idx <- grep(paste0("^", human_prefix), rownames(counts))
  m_idx <- grep(paste0("^", mouse_prefix), rownames(counts))
  if (length(h_idx) == 0) {
    stop("No ", human_prefix, " features in ", unique(obj$orig.ident), call. = FALSE)
  }
  h <- Matrix::colSums(counts[h_idx, , drop = FALSE])
  m <- if (length(m_idx) > 0) Matrix::colSums(counts[m_idx, , drop = FALSE]) else rep(0, ncol(counts))
  frac <- h / pmax(h + m, 1)
  ifelse(frac >= 0.9, "GRCh38", "not_human")
}

add_percent_mt <- function(obj) {
  mt_genes <- grep(paste0("^", human_prefix, "MT-"), rownames(obj), value = TRUE)
  if (length(mt_genes) == 0) {
    log_msg("  warning: no ", human_prefix, "MT- genes; percent.mt = 0")
    obj$percent.mt <- 0
  } else {
    obj$percent.mt <- PercentageFeatureSet(obj, features = mt_genes)
  }
  obj
}

run_doubletfinder <- function(obj, label) {
  n <- n_cells(obj)
  if (n < 100) {
    log_msg("  ", label, ": <100 cells, skip DoubletFinder (all Singlet)")
    obj$DF_Classification <- "Singlet"
    obj$DF_pANN <- 0
    return(obj)
  }
  # Literal F06_01 setting. DoubletFinder rebuilds the artificial-doublet
  # object with SCTransform and RunPCA(npcs=length(PCs)), so PCs 1:75 do
  # not depend on the 30-PC reduction used above for clustering.
  pcs <- seq_len(min(75L, n - 1L, n_features(obj) - 1L))
  if (length(pcs) < 5) stop(label, " has <5 PCs for DoubletFinder", call. = FALSE)

  # F06_01: estimated_doublet_rate=0.076, pN=0.25, paramSweep_v3/doubletFinder_v3.
  # Current DoubletFinder dropped the _v3 suffix; all statistical settings
  # remain as published, including sct=TRUE.
  sweep.res <- paramSweep(obj, PCs = pcs, sct = TRUE, num.cores = n_threads)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  if (!is.finite(pK) || pK <= 0) {
    stop(label, ": published DoubletFinder pK sweep produced no valid pK",
         call. = FALSE)
  }
  log_msg(sprintf("  %s DoubletFinder pK=%.3f (published sweep)", label, pK))

  annotations <- obj$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.076 * n)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  obj <- doubletFinder(
    obj, PCs = pcs, pN = 0.25, pK = pK, nExp = nExp_poi.adj,
    reuse.pANN = NULL, sct = TRUE
  )
  md <- obj[[]]
  df_cols <- grep("^DF\\.classifications", colnames(md), value = TRUE)
  pann_cols <- grep("^pANN", colnames(md), value = TRUE)
  if (length(df_cols) == 0) stop(label, ": DoubletFinder did not add DF.classifications", call. = FALSE)
  obj$DF_Classification <- md[[df_cols[1]]]
  obj$DF_pANN <- if (length(pann_cols) > 0) md[[pann_cols[1]]] else NA_real_
  log_msg(sprintf("  %s: %d singlets, %d doublets", label,
                  sum(obj$DF_Classification == "Singlet"),
                  sum(obj$DF_Classification == "Doublet")))
  obj
}

write_qc_table <- function(obj_list, file) {
  dt <- rbindlist(lapply(names(obj_list), function(nm) {
    md <- obj_list[[nm]][[]]
    data.table(
      sample = nm,
      nFeature_RNA = md$nFeature_RNA,
      nCount_RNA = md$nCount_RNA,
      percent.mt = md$percent.mt
    )
  }))
  fwrite(dt, file, sep = "\t")
}

# ============================================================================
log_msg("Seurat ", as.character(packageVersion("Seurat")),
        " / SeuratObject ", as.character(packageVersion("SeuratObject")),
        " / DoubletFinder ", as.character(packageVersion("DoubletFinder")))
log_msg("threads=", n_threads, "  FORCE=", force_rerun)
preflight_cellranger()


# ---------------------------- Phase D ---------------------------------------
if (!force_rerun && file.exists(phase_d_rda)) {
  log_msg("\nPhase D: resume from ", phase_d_rda)
  load(phase_d_rda)
  if (!exists("KO") || length(KO) != 4) stop("Corrupt Phase D checkpoint", call. = FALSE)
} else {
  log_msg("\n============================================================")
  log_msg("Phase D: Seurat preprocessing of ADAR-KO samples")
  log_msg("============================================================")

  ter_KO.SO.list <- list()
  for (i in seq_along(samples)) {
    raw <- Read10X(file.path(cr_dir, samples[i], "outs", "filtered_feature_bc_matrix"))
    mat <- extract_counts(raw)
    rownames(mat) <- gsub("_", "-", rownames(mat), fixed = TRUE)
    log_msg(sprintf("  %s Read10X: %d genes x %d cells", samples[i],
                    n_features(mat), n_cells(mat)))
    if (ncol(mat) == 0) stop(samples[i], " count matrix has 0 cells", call. = FALSE)

    obj <- new_seurat(mat, project = "KO_H1")
    # F06_01: min.cells = round(0.001 * ncol); min.features = 50 then subset >200
    obj <- filter_min_cells(obj, 0.001)
    keep50 <- WhichCells(obj, expression = nFeature_RNA >= 50)
    if (length(keep50) < 200) {
      stop(ter_list[i], " after min.features=50: ", length(keep50), " cells",
           call. = FALSE)
    }
    obj <- subset(obj, cells = keep50)
    obj$teratoma <- ter_list[i]
    obj <- RenameCells(obj, new.names = paste0(ter_list[i], "_", colnames(obj)))
    obj <- add_percent_mt(obj)
    log_msg(sprintf("  %s loaded: %d cells  median nFeature=%.0f  mean percent.mt=%.2f",
                    ter_list[i], n_cells(obj), median(obj$nFeature_RNA), mean(obj$percent.mt)))
    ter_KO.SO.list[[i]] <- obj
  }
  names(ter_KO.SO.list) <- ter_list
  write_qc_table(ter_KO.SO.list, file.path(out_dir, "ko_qc_prefilter.tsv"))

  for (i in seq_along(ter_KO.SO.list)) {
    n_before <- n_cells(ter_KO.SO.list[[i]])
    keep <- WhichCells(ter_KO.SO.list[[i]],
                       expression = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
    if (length(keep) < 200) {
      stop(ter_list[i], " QC left ", length(keep), " / ", n_before,
           " cells — check percent.mt pattern and CellRanger quality", call. = FALSE)
    }
    ter_KO.SO.list[[i]] <- subset(ter_KO.SO.list[[i]], cells = keep)
    log_msg(sprintf("  %s QC: %d → %d cells", ter_list[i], n_before,
                    n_cells(ter_KO.SO.list[[i]])))
  }

  for (i in seq_along(ter_KO.SO.list)) {
    obj <- ter_KO.SO.list[[i]]
    n_pcs <- min(30L, n_cells(obj) - 1L, n_features(obj) - 1L)
    # F06_01: LogNormalize 1e4, vst 4000, PCA 30, UMAP 1:30, k=10, res=0.6, algorithm=1
    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = min(4000L, nrow(obj) - 1L), verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:n_pcs, verbose = FALSE)
    obj <- FindNeighbors(obj, k.param = 10, reduction = "pca",
                         dims = 1:n_pcs, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.6, algorithm = 1, verbose = FALSE)
    ter_KO.SO.list[[i]] <- run_doubletfinder(obj, ter_list[i])
  }

  for (i in seq_along(ter_KO.SO.list)) {
    obj <- ter_KO.SO.list[[i]]
    gem_calls <- read_gem_calls(samples[i], ter_list[i], colnames(obj))
    if (!is.null(gem_calls) && length(gem_calls) >= 0.5 * ncol(obj)) {
      obj <- obj[, names(gem_calls)]
      obj$cr_call <- gem_calls[colnames(obj)]
      log_msg(sprintf("  %s species from gem_classification: %s",
                      ter_list[i], paste(names(table(obj$cr_call)), table(obj$cr_call),
                                         sep = "=", collapse = ", ")))
    } else {
      obj$cr_call <- human_frac_calls(obj)
      log_msg(sprintf("  %s species from human UMI fraction >= 0.9: %s",
                      ter_list[i], paste(names(table(obj$cr_call)), table(obj$cr_call),
                                         sep = "=", collapse = ", ")))
    }

    keep <- colnames(obj)[obj$DF_Classification == "Singlet" & obj$cr_call == "GRCh38"]
    if (length(keep) < 200) {
      stop(ter_list[i], " after singlet + human filter: ", length(keep), " cells",
           call. = FALSE)
    }
    ter_KO.SO.list[[i]] <- subset(obj, cells = keep)
    log_msg(sprintf("  %s after singlet + human filter: %d cells",
                    ter_list[i], n_cells(ter_KO.SO.list[[i]])))
  }

  KO <- list()
  for (i in seq_along(ter_KO.SO.list)) {
    obj <- ter_KO.SO.list[[i]]
    genes.human <- grep(paste0("^", human_prefix), rownames(obj), value = TRUE)
    if (length(genes.human) < 500) {
      stop(ter_list[i], " has only ", length(genes.human), " human genes", call. = FALSE)
    }
    counts.temp <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes.human, , drop = FALSE]
    rownames(counts.temp) <- sub(paste0("^", human_prefix), "", rownames(counts.temp))
    genes.keep <- rownames(counts.temp)[!grepl("^RP11|^MT-", rownames(counts.temp))]
    counts.clean <- counts.temp[genes.keep, , drop = FALSE]
    KO[[i]] <- new_seurat(counts.clean, meta = obj[[]], project = "KO_H1")
    require_ncells(KO[[i]], paste(ter_list[i], "human-only object"), 200)
  }
  names(KO) <- ter_list
  save(KO, file = phase_d_rda)
  log_msg(sprintf("Phase D complete. Total human singlet cells: %d",
                  sum(vapply(KO, n_cells, 1L))))
}


# ---------------------------- Phase E ---------------------------------------
if (!force_rerun && file.exists(phase_e_rds)) {
  log_msg("\nPhase E: resume from ", phase_e_rds)
} else {
  log_msg("\n============================================================")
  log_msg("Phase E: Build WT teratoma reference")
  log_msg("============================================================")

  counts <- Read10X(file.path(wt_dir, "dm-ter-human-merged"))
  doublets_file <- file.path(wt_dir, "GSE156170_dm-ter-merged_doublets.txt")
  require_file(doublets_file, "WT doublet flags")
  doublets <- as.numeric(scan(doublets_file, sep = "\n", what = numeric(), quiet = TRUE))
  doublets[is.na(doublets)] <- 1
  if (length(doublets) != ncol(counts)) {
    stop("doublet flags (", length(doublets), ") != WT cells (", n_cells(counts), ")",
         call. = FALSE)
  }
  log_msg(sprintf("  Doublet rate: %.1f%% (%d/%d)",
                  100 * sum(doublets > 0) / length(doublets),
                  sum(doublets > 0), length(doublets)))
  counts <- counts[, doublets == 0]

  metadata <- fread(file.path(wt_dir, "dm-ter-human-merged_metadata.tsv"))
  metadata[, cell_key := gsub("\\.", "-", cell)]
  colnames(counts) <- gsub("\\.", "-", colnames(counts))
  keep_cells <- intersect(colnames(counts), metadata$cell_key)
  log_msg(sprintf("  Cells with published labels after doublet filter: %d / %d",
                  length(keep_cells), n_cells(counts)))
  if (length(keep_cells) < 1000) {
    stop("Too few WT cells matched metadata (barcode '-' vs '.' mismatch?)", call. = FALSE)
  }
  counts <- counts[, keep_cells]
  setkey(metadata, cell_key)

  rownames(counts) <- gsub("^hg19_", "", rownames(counts))
  rownames(counts) <- gsub("^mm10_", "", rownames(counts))
  # Literal yanwu Figure1/01_human_clustering.R filter.
  genes.keep <- rownames(counts)[!grepl("RP11|MT", rownames(counts))]
  counts <- counts[genes.keep, ]

  suffix <- sub(".*[-.]", "", colnames(counts))
  batch <- factor(paste0("ter", suffix),
                  levels = paste0("ter", sort(unique(as.integer(suffix)))))
  names(batch) <- colnames(counts)
  barcoded <- ifelse(as.character(batch) %in% paste0("ter", 1:4), "nonbc", "bc")
  names(barcoded) <- colnames(counts)
  meta.data <- data.frame(batch = batch, barcoded = factor(barcoded),
                          row.names = colnames(counts))

  obj <- new_seurat(counts, meta = meta.data, project = "WT")
  # yanwu 01_human_clustering.R: min.features=200, min.cells=round(0.001*ncol)
  keep <- WhichCells(obj, expression = nFeature_RNA >= 200)
  obj <- subset(obj, cells = keep)
  obj <- filter_min_cells(obj, 0.001)
  require_ncells(obj, "WT after min.features=200", 1000)
  log_msg(sprintf("  WT object: %d cells x %d genes across %d batches",
                  n_cells(obj), n_features(obj), length(unique(obj$batch))))

  obj.list <- SplitObject(obj, split.by = "batch")
  rm(obj, counts)
  invisible(gc())

  log_msg("  Normalizing and finding variable features per batch...")
  for (i in seq_along(obj.list)) {
    obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
    obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst",
                                          nfeatures = 4000, verbose = FALSE)
  }

  log_msg("  Finding integration anchors...")
  obj.anchors <- FindIntegrationAnchors(object.list = obj.list,
                                        anchor.features = 4000, dims = 1:30)
  log_msg("  Integrating data...")
  obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)
  DefaultAssay(obj.integrated) <- "integrated"
  # yanwu ScaleData on all integrated features; return.model for F06_01 MapQuery
  obj.integrated <- ScaleData(obj.integrated,
                              features = rownames(GetAssayData(obj.integrated)),
                              verbose = FALSE)
  obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:30,
                            metric = "cosine", return.model = TRUE, verbose = FALSE)
  obj.integrated <- FindNeighbors(obj.integrated, k.param = 10,
                                  reduction = "pca", dims = 1:30)
  obj.integrated <- FindClusters(obj.integrated, resolution = 0.4,
                                 algorithm = 2, verbose = FALSE)

  cl <- metadata[colnames(obj.integrated), cluster]
  obj.integrated$cluster <- cl
  labeled <- colnames(obj.integrated)[!is.na(obj.integrated$cluster)]
  if (length(labeled) < n_cells(obj.integrated)) {
    log_msg(sprintf("  Dropping %d unlabeled WT cells before transfer",
                    n_cells(obj.integrated) - length(labeled)))
    obj.integrated <- subset(obj.integrated, cells = labeled)
  }
  Idents(obj.integrated) <- obj.integrated$cluster
  log_msg("  WT reference cell type distribution:")
  print(sort(table(obj.integrated$cluster), decreasing = TRUE))

  umap <- Embeddings(obj.integrated, "umap")
  fwrite(data.table(
    cell = rownames(umap),
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    cluster = as.character(obj.integrated$cluster),
    batch = as.character(obj.integrated$batch)
  ), file.path(out_dir, "wt_umap.tsv"), sep = "\t")

  saveRDS(obj.integrated, phase_e_rds)
  rm(obj.list, obj.anchors, obj.integrated)
  invisible(gc())
  log_msg("Phase E complete. WT reference saved.")
}


# ---------------------------- Phase F ---------------------------------------
if (!force_rerun && file.exists(phase_f_rda)) {
  log_msg("\nPhase F: resume from ", phase_f_rda)
  load(phase_f_rda)
  if (!exists("KO.query") || length(KO.query) != 4) {
    stop("Corrupt Phase F checkpoint", call. = FALSE)
  }
} else {
  log_msg("\n============================================================")
  log_msg("Phase F: Label transfer WT → ADAR-KO")
  log_msg("============================================================")

  wt <- readRDS(phase_e_rds)
  if (is.null(wt$cluster) || all(is.na(wt$cluster))) {
    stop("WT reference has no cluster labels — Phase E barcode match failed", call. = FALSE)
  }
  if (!"umap" %in% Reductions(wt)) stop("WT reference missing UMAP model", call. = FALSE)

  KO.query <- list()
  for (i in seq_along(KO)) {
    lab <- names(KO)[i]
    log_msg("  Transferring labels for ", lab, "...")
    # F06_01 label transfer (FindTransferAnchors + TransferData + MapQuery)
    q <- NormalizeData(KO[[i]], verbose = FALSE)
    anchors <- FindTransferAnchors(
      reference = wt, query = q,
      dims = 1:30, reference.reduction = "pca"
    )
    predictions <- TransferData(
      anchorset = anchors, refdata = wt$cluster, dims = 1:30
    )
    q <- AddMetaData(q, metadata = predictions)
    q <- NormalizeData(q, verbose = FALSE)
    q <- FindVariableFeatures(q, selection.method = "vst", nfeatures = 2000,
                              verbose = FALSE)
    configure_seurat_future()
    q <- MapQuery(
      anchorset = anchors, reference = wt, query = q,
      refdata = list(celltype = "cluster"),
      reference.reduction = "pca", reduction.model = "umap", verbose = FALSE
    )
    if (is.null(q$predicted.celltype) && "predicted.id" %in% colnames(q[[]])) {
      q$predicted.celltype <- as.character(q$predicted.id)
    }
    if (is.null(q$predicted.celltype)) {
      stop(lab, ": MapQuery/TransferData did not add predicted.celltype",
           call. = FALSE)
    }
    log_msg(sprintf("  %s (%d cells); min score=%.3f",
                    lab, n_cells(q), min(q$prediction.score.max, na.rm = TRUE)))
    print(sort(table(q$predicted.celltype), decreasing = TRUE))
    KO.query[[i]] <- q
  }
  names(KO.query) <- names(KO)
  save(KO.query, file = phase_f_rda)
  rm(wt)
  invisible(gc())
  log_msg("Phase F complete.")
}


# ---------------------------- Phase G ---------------------------------------
log_msg("\n============================================================")
log_msg("Phase G: Export barcode → celltype TSVs")
log_msg("============================================================")

sample_map <- c(ter1 = "T1", ter2 = "T2", ter3 = "T3", ter4 = "T4")
ss_rows <- list()

for (i in seq_along(ter_list)) {
  obj <- KO.query[[i]]
  sample_id <- sample_map[[ter_list[i]]]
  raw_barcodes <- colnames(obj)
  clean_barcodes <- sub("^ter[0-9]+_", "", raw_barcodes)
  clean_barcodes <- sub("-[0-9]+$", "", clean_barcodes)

  df <- data.table(
    barcode = clean_barcodes,
    celltype = sanitize_celltype(obj$predicted.celltype),
    score = obj$prediction.score.max
  )
  if (anyDuplicated(df$barcode)) {
    log_msg("  warning: ", sample_id, " has duplicated barcodes after stripping GEM suffix")
  }

  ct_file <- file.path(pm_dir, paste0("barcodes_celltype_", sample_id, ".tsv"))
  fwrite(df[, .(barcode, celltype)], ct_file, sep = "\t", quote = FALSE, col.names = FALSE)

  scored_file <- file.path(out_dir, paste0("barcodes_celltype_scored_", sample_id, ".tsv"))
  fwrite(df, scored_file, sep = "\t", quote = FALSE)

  bc_file <- file.path(pm_dir, paste0("barcodes_qcfiltered_", sample_id, ".tsv"))
  fwrite(df[, .(barcode)], bc_file, sep = "\t", quote = FALSE, col.names = FALSE)

  log_msg(sprintf("  %s: %d cells → %s (%d types)",
                  sample_id, nrow(df), basename(ct_file), uniqueN(df$celltype)))

  bam <- file.path(pm_dir, "01_starsolobams", paste0(sample_id, "_merged.bam"))
  ss_rows[[i]] <- data.table(
    sample_id = sample_id,
    bam = bam,
    celltype_annot = ct_file,
    called_barcodes = bc_file,
    platform = "10xv3",
    fastq_r1 = "",
    fastq_r2 = ""
  )
}

pred_plot <- rbindlist(lapply(seq_along(ter_list), function(i) {
  obj <- KO.query[[i]]
  data.table(
    sample = sample_map[[ter_list[i]]],
    celltype = sanitize_celltype(obj$predicted.celltype),
    score = obj$prediction.score.max
  )
}))
fwrite(pred_plot, file.path(out_dir, "ko_predicted_celltypes.tsv"), sep = "\t")

umap_plot <- rbindlist(lapply(seq_along(ter_list), function(i) {
  obj <- KO.query[[i]]
  red <- intersect(c("ref.umap", "umap"), Reductions(obj))[1]
  if (is.na(red)) return(NULL)
  emb <- Embeddings(obj, red)
  data.table(
    sample = sample_map[[ter_list[i]]],
    cell = colnames(obj),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    celltype = sanitize_celltype(obj$predicted.celltype)
  )
}), fill = TRUE)
if (nrow(umap_plot) > 0L) {
  fwrite(umap_plot, file.path(out_dir, "ko_query_umap.tsv"), sep = "\t")
}

ss <- rbindlist(ss_rows)
ss_file <- file.path(pm_dir, "samplesheet_pm.csv")
fwrite(ss, ss_file)
log_msg("Wrote ", ss_file)
if (any(!file.exists(ss$bam))) {
  log_msg("  note: merged BAM(s) not on disk yet; samplesheet paths are ready")
}

log_msg("\nAll phases complete.")
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
