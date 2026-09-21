#!/usr/bin/env Rscript
# goodCells = annotated ∪ QC-passing unannotated
#
# Candidate pool = STARsolo/CellRanger filtered barcodes ∪ annotated barcodes
# (or called_barcodes ∪ annotated when no count matrix is available).
# Annotated cells are assumed QC-passed. Only unannotated cells are QC-filtered.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args = function(args) {
  out = list(
    sample = "",
    annot = "",
    called = "",
    qc_barcodes = "",
    sample_cells = "",
    matrix_dir = "",
    gtf = "",
    mito_chr = "",
    qc_min_features = 200L,
    qc_max_features = 6000L,
    qc_max_pct_mt = 20,
    qc_min_counts = 500L,
    output_barcodes = "goodcells_barcodes.tsv",
    output_summary = "barcode_qc_summary.tsv",
    output_metrics = "unannotated_qc_metrics.tsv"
  )
  i = 1L
  while (i <= length(args)) {
    key = sub("^--", "", args[[i]])
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for --", key)
    }
    val = args[[i + 1L]]
    if (key %in% names(out)) {
      if (is.integer(out[[key]])) {
        out[[key]] = as.integer(val)
      } else if (is.numeric(out[[key]])) {
        out[[key]] = as.numeric(val)
      } else {
        out[[key]] = val
      }
    }
    i = i + 2L
  }
  out
}

opt = parse_args(commandArgs(trailingOnly = TRUE))

if (!nzchar(opt$sample)) stop("--sample is required")
if (!nzchar(opt$annot) || !file.exists(opt$annot)) stop("--annot file is required")
if (!nzchar(opt$gtf) || !file.exists(opt$gtf)) stop("--gtf file is required")

bc_key = function(x) {
  x = trimws(as.character(x))
  x = x[nzchar(x)]
  sub("-\\d+$", "", x)
}

read_barcode_col1 = function(path) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(character())
  }
  dt = fread(path, header = FALSE, sep = "\t", select = 1L, col.names = "barcode")
  unique(trimws(as.character(dt$barcode)))
}

display_for_keys = function(keys, preferred, fallback) {
  pref_dt = data.table(barcode = unique(preferred), keyv = bc_key(unique(preferred)))
  fall_dt = data.table(barcode = unique(fallback), keyv = bc_key(unique(fallback)))
  pref_dt = pref_dt[!duplicated(keyv)]
  fall_dt = fall_dt[!duplicated(keyv)]
  out = pref_dt[data.table(keyv = keys), on = "keyv"]
  miss = is.na(out$barcode)
  if (any(miss)) {
    fill = fall_dt[data.table(keyv = keys[miss]), on = "keyv"]
    out$barcode[miss] = fill$barcode
  }
  out$barcode[is.na(out$barcode)] = keys[is.na(out$barcode)]
  out$barcode
}

DEFAULT_MITO_CHRS = c(
  "chrM", "chrMT", "MT", "M", "mitochondrion_genome",
  "chrMito", "mito", "mitochondrion"
)

detect_mito_chrs = function(gtf, override) {
  if (nzchar(override)) {
    return(unique(trimws(unlist(strsplit(override, ",", fixed = TRUE)))))
  }
  cmd = sprintf(
    "awk -F'\\t' 'BEGIN{IGNORECASE=1} $1 ~ /^(chrM|chrMT|MT|M|mitochondrion_genome|chrMito|mito)$/ || $1 ~ /mito/ {print $1}' %s | sort -u",
    shQuote(gtf)
  )
  found = unique(trimws(system(cmd, intern = TRUE)))
  found = found[nzchar(found)]
  if (length(found) == 0) {
    warning("No mitochondrial chromosomes detected in GTF; pct_mt will be 0. Use --mito_chr to override.")
  }
  found
}

read_mito_genes = function(gtf, mito_chrs) {
  if (length(mito_chrs) == 0) return(character())
  cond = paste(sprintf("$1==\"%s\"", mito_chrs), collapse = " || ")
  cmd = sprintf("awk -F'\\t' '%s {print $9}' %s", cond, shQuote(gtf))
  attrs = system(cmd, intern = TRUE)
  genes = character()
  if (length(attrs) > 0) {
    gid = regmatches(attrs, regexpr('gene_id "[^"]+"', attrs))
    gnm = regmatches(attrs, regexpr('gene_name "[^"]+"', attrs))
    genes = c(
      sub('gene_id "', "", sub('"$', "", gid)),
      sub('gene_name "', "", sub('"$', "", gnm))
    )
  }
  unique(genes[nzchar(genes)])
}

read_one_matrix = function(matrix_dir) {
  mtx = file.path(matrix_dir, "matrix.mtx")
  if (!file.exists(mtx)) mtx = file.path(matrix_dir, "matrix.mtx.gz")
  feat = file.path(matrix_dir, "features.tsv")
  if (!file.exists(feat)) feat = file.path(matrix_dir, "features.tsv.gz")
  if (!file.exists(feat)) feat = file.path(matrix_dir, "genes.tsv")
  if (!file.exists(feat)) feat = file.path(matrix_dir, "genes.tsv.gz")
  bc = file.path(matrix_dir, "barcodes.tsv")
  if (!file.exists(bc)) bc = file.path(matrix_dir, "barcodes.tsv.gz")
  if (!file.exists(mtx) || !file.exists(feat) || !file.exists(bc)) {
    stop("Count matrix directory is missing matrix/features/barcodes: ", matrix_dir)
  }

  con = file(mtx, "r")
  on.exit(close(con), add = TRUE)
  skip = 0L
  dim_line = ""
  repeat {
    line = readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("Empty MTX file: ", mtx)
    if (!startsWith(line, "%")) {
      dim_line = line
      break
    }
    skip = skip + 1L
  }
  dims = as.integer(strsplit(dim_line, "\\s+")[[1]])
  if (length(dims) < 3L) stop("Invalid MTX header: ", mtx)

  counts = fread(mtx, skip = skip + 1L, header = FALSE,
                 col.names = c("gene_i", "cell_i", "count"))
  features = fread(feat, header = FALSE, sep = "\t")
  barcodes = fread(bc, header = FALSE, sep = "\t", select = 1L, col.names = "barcode")
  gene_id = trimws(as.character(features[[1]]))
  gene_name = if (ncol(features) >= 2L) trimws(as.character(features[[2]])) else gene_id
  barcodes[, barcode := trimws(barcode)]

  counts[, gene_id := gene_id[gene_i]]
  counts[, gene_name := gene_name[gene_i]]
  counts[, barcode := barcodes$barcode[cell_i]]
  counts[, c("gene_i", "cell_i") := NULL]
  counts
}

load_matrices = function(matrix_dirs) {
  dirs = trimws(unlist(strsplit(matrix_dirs, ",", fixed = TRUE)))
  dirs = dirs[nzchar(dirs)]
  if (length(dirs) == 0L) return(NULL)
  parts = lapply(dirs, read_one_matrix)
  dt = rbindlist(parts, use.names = TRUE, fill = TRUE)
  dt[, .(count = sum(count)), by = .(gene_id, gene_name, barcode)]
}

annot_bc = read_barcode_col1(opt$annot)
if (length(annot_bc) == 0L) stop("No annotated barcodes in ", opt$annot)

called_bc = read_barcode_col1(opt$called)
qc_override = read_barcode_col1(opt$qc_barcodes)
sample_bc = read_barcode_col1(opt$sample_cells)

if (length(sample_bc) > 0L) {
  sample_keys = unique(bc_key(sample_bc))
  annot_bc = annot_bc[bc_key(annot_bc) %in% sample_keys]
  called_bc = called_bc[bc_key(called_bc) %in% sample_keys]
  if (length(annot_bc) == 0L) {
    warning("No annotated barcodes overlap sample cells for ", opt$sample)
  }
}

annot_keys = unique(bc_key(annot_bc))
sample_keys = unique(bc_key(sample_bc))

mito_chrs = detect_mito_chrs(opt$gtf, opt$mito_chr)
mito_genes = read_mito_genes(opt$gtf, mito_chrs)

qc_source = ""
metrics = data.table()
unannotated_keys = character()
good_unannot_keys = character()
fail_unannot_keys = character()
matrix_keys = character()
n_matrix_barcodes = 0L

if (length(qc_override) > 0L) {
  qc_source = "qc_barcodes"
  good_keys = unique(c(annot_keys, bc_key(qc_override)))
} else {
  counts = load_matrices(opt$matrix_dir)
  if (!is.null(counts)) {
    qc_source = "matrix"
    matrix_bc = unique(counts$barcode)
    n_matrix_barcodes = length(matrix_bc)
    matrix_keys = unique(bc_key(matrix_bc))
    pool_keys = unique(c(matrix_keys, annot_keys))
    unannotated_keys = setdiff(pool_keys, annot_keys)

    cell_metrics = counts[, .(
      nFeature = uniqueN(gene_id[count > 0]),
      nCount = sum(count),
      nCount_mt = sum(count[gene_id %in% mito_genes | gene_name %in% mito_genes])
    ), by = barcode]
    cell_metrics[, keyv := bc_key(barcode)]
    cell_metrics[, pct_mt := fifelse(nCount > 0, 100 * nCount_mt / nCount, 0)]
    # Keep one row per normalized barcode (sum already done at count level)
    cell_metrics = cell_metrics[, .(
      barcode = barcode[1],
      nFeature = max(nFeature),
      nCount = sum(nCount),
      nCount_mt = sum(nCount_mt),
      pct_mt = {
        tot = sum(nCount)
        fifelse(tot > 0, 100 * sum(nCount_mt) / tot, 0)
      }
    ), by = keyv]

    un_metrics = cell_metrics[keyv %in% unannotated_keys]
    # Barcodes present only in annotation (not in matrix) have no metrics;
    # they are treated as annotated and already excluded from unannotated_keys.
    un_metrics[, qc_pass :=
      nFeature >= opt$qc_min_features &
      nFeature <= opt$qc_max_features &
      nCount >= opt$qc_min_counts &
      pct_mt <= opt$qc_max_pct_mt
    ]
    # Unannotated barcodes missing from the matrix cannot be QC'd — drop them
    missing_un = setdiff(unannotated_keys, un_metrics$keyv)
    good_unannot_keys = un_metrics[qc_pass == TRUE, keyv]
    fail_unannot_keys = unique(c(un_metrics[qc_pass == FALSE, keyv], missing_un))
    metrics = un_metrics
    good_keys = unique(c(annot_keys, good_unannot_keys))
  } else if (length(sample_keys) > 0L || length(called_bc) > 0L) {
    qc_source = if (length(sample_keys) > 0L) "sample_cells" else "called_fallback"
    pool_keys = unique(c(sample_keys, bc_key(called_bc), annot_keys))
    unannotated_keys = setdiff(pool_keys, annot_keys)
    # No matrix: keep caller/samplesheet cells (Smart-seq2 QC already applied)
    good_unannot_keys = unannotated_keys
    fail_unannot_keys = character()
    good_keys = pool_keys
  } else {
    qc_source = "annotated_only"
    warning("No count matrix, called_barcodes, or qc_barcodes for ", opt$sample,
            "; goodCells will equal annotated.")
    good_keys = annot_keys
  }
}

good_display = display_for_keys(good_keys, annot_bc, c(called_bc, qc_override, metrics$barcode))
good_dt = data.table(barcode = good_display, label = "goodCells")
fwrite(good_dt, opt$output_barcodes, sep = "\t", col.names = FALSE)

if (nrow(metrics) > 0L) {
  metrics_out = metrics[, .(
    sample = opt$sample,
    barcode,
    nFeature, nCount, nCount_mt, pct_mt, qc_pass
  )]
} else {
  metrics_out = data.table(
    sample = character(), barcode = character(),
    nFeature = integer(), nCount = numeric(), nCount_mt = numeric(),
    pct_mt = numeric(), qc_pass = logical()
  )
}
fwrite(metrics_out, opt$output_metrics, sep = "\t")

summary_dt = data.table(
  sample = opt$sample,
  qc_source = qc_source,
  n_annotated = length(annot_keys),
  n_matrix_barcodes = n_matrix_barcodes,
  n_unannotated = length(unannotated_keys),
  n_unannotated_qc_pass = length(good_unannot_keys),
  n_unannotated_qc_fail = length(fail_unannot_keys),
  n_goodcells = length(good_keys),
  n_mito_genes = length(mito_genes),
  mito_chroms = paste(mito_chrs, collapse = ","),
  qc_min_features = opt$qc_min_features,
  qc_max_features = opt$qc_max_features,
  qc_min_counts = opt$qc_min_counts,
  qc_max_pct_mt = opt$qc_max_pct_mt
)
fwrite(summary_dt, opt$output_summary, sep = "\t")

cat(sprintf(
  "[%s] source=%s annotated=%d unannotated=%d pass=%d fail=%d goodCells=%d mito_genes=%d mito_chr=%s\n",
  opt$sample, qc_source, length(annot_keys), length(unannotated_keys),
  length(good_unannot_keys), length(fail_unannot_keys), length(good_keys),
  length(mito_genes), paste(mito_chrs, collapse = ",")
))
