#!/usr/bin/env Rscript
# Validate a pipeline samplesheet against nextflow/SAMPLESHEET_SPEC.md.
#
#   Rscript validate_samplesheet.R <samplesheet.csv> [--min-cells 50] [--quiet]
#
# Exit 0 = conforms (warnings allowed), 1 = one or more errors.
# Reports the group/sample/celltype structure the cell-type split derives, so a
# dataset can be checked before any splitting is launched.

suppressPackageStartupMessages(library(data.table))

# ── args ─────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] %in% c("-h", "--help")) {
    cat("Usage: validate_samplesheet.R <samplesheet.csv>",
        "[--min-cells N] [--dataset NAME] [--tsv-out PATH] [--emit-plan DIR] [--quiet]\n")
    cat("  --min-cells  minimum cells per cell type. DEFAULT 0 = keep all.\n")
    cat("               Any non-zero value is a filter not in the original analysis\n")
    cat("               and must be recorded in the scope TSV before use.\n")
    cat("  --emit-plan  write the machine-readable cell-type split plan to DIR\n")
    cat("               (_expected.txt, _plan_<unit>.tsv, _cells_<unit>.tsv,\n")
    cat("               _keep_<unit>.txt). Written only if the sheet validates.\n")
    quit(status = if (length(args)) 0L else 1L)
}
sheet_path <- args[1]
getopt <- function(flag, default) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) default else args[i + 1L]
}
# DEFAULT 0 — no filter. A cell-type minimum is not part of the original
# analysis, so it is opt-in and recorded, never silently applied.
MIN_CELLS <- as.integer(getopt("--min-cells", "0"))
DATASET   <- getopt("--dataset", basename(dirname(normalizePath(sheet_path, mustWork = FALSE))))
TSV_OUT   <- getopt("--tsv-out", NA_character_)
# --emit-plan writes the machine-readable split plan (see "split plan" below).
# It is what makes split_gt_bams_celltype.sh dataset-agnostic: the splitter reads
# these files instead of re-deriving units, cell types or cell counts, so it
# cannot drift from what was validated and recorded in CELLTYPE_SCOPE.tsv.
PLAN_DIR  <- getopt("--emit-plan", NA_character_)
QUIET     <- "--quiet" %in% args
scope_rows <- list()
plan_data  <- list()

PLATFORMS_BARCODED <- c("10xv2", "10xv3", "celseq2")
PLATFORMS_PLATE    <- c("smartseq2")
PLATFORMS_ALL      <- c(PLATFORMS_BARCODED, PLATFORMS_PLATE)
REQUIRED           <- c("sample_id", "bam", "celltype_annot", "called_barcodes", "platform")

errs <- character(0)
warns <- character(0)
err  <- function(...) errs  <<- c(errs,  paste0(...))
warn <- function(...) warns <<- c(warns, paste0(...))
say  <- function(...) if (!QUIET) cat(..., "\n", sep = "")

# ── load ─────────────────────────────────────────────────────────────────────
if (!file.exists(sheet_path)) {
    cat("ERROR: no such samplesheet: ", sheet_path, "\n", sep = "")
    quit(status = 1L)
}
ss <- tryCatch(
    fread(sheet_path, colClasses = "character", na.strings = c("", "NA")),
    error = function(e) { cat("ERROR: cannot parse CSV: ", conditionMessage(e), "\n", sep = ""); quit(status = 1L) }
)
say("samplesheet : ", sheet_path)
say("rows        : ", nrow(ss))

# 0. duplicated column names. fread keeps BOTH, and `ss$group` then silently
# resolves to the FIRST — so a sheet with `group` twice validated clean while
# being malformed (exactly what happened to samplesheet_pm_gt_persample.csv,
# where the group-adding edit ran twice). Checked before anything else, because
# every later column lookup is ambiguous until it is ruled out.
dup_cols <- unique(names(ss)[duplicated(names(ss))])
if (length(dup_cols)) {
    err("duplicated column name(s): ", paste(dup_cols, collapse = ", "),
        " — every value lookup silently takes the first, so the sheet must be fixed, not tolerated")
    cat("\n", paste0("ERROR: ", errs, collapse = "\n"), "\n", sep = "")
    quit(status = 1L)
}

# 1. required columns
missing_cols <- setdiff(REQUIRED, names(ss))
if (length(missing_cols)) {
    err("missing required column(s): ", paste(missing_cols, collapse = ", "))
    cat("\n", paste0("ERROR: ", errs, collapse = "\n"), "\n", sep = "")
    quit(status = 1L)
}

# Optional columns normalised so downstream logic is uniform.
if (!"group"   %in% names(ss)) ss[, group   := NA_character_]
if (!"cell_id" %in% names(ss)) ss[, cell_id := NA_character_]
ss[is.na(group) | !nzchar(trimws(group)), group := "all"]
ss[, group := trimws(group)]

# `group` and `sample_id` become path components in the split plan and in the
# emitted BAM names, so they must be filename-safe. Rejected rather than
# sanitised: silently renaming a pooling unit would make the scope table
# disagree with the sample_id the pipeline then reports.
for (cl in c("group", "sample_id")) {
    bad <- unique(ss[[cl]][!grepl("^[A-Za-z0-9._-]+$", ss[[cl]])])
    if (length(bad)) err(length(bad), " ", cl, " value(s) contain characters outside ",
                         "[A-Za-z0-9._-] and cannot be used in a filename: ",
                         paste(head(bad, 3), collapse = ", "))
}

# 2. no blanks in required fields
for (cl in REQUIRED) {
    n <- ss[is.na(get(cl)) | !nzchar(trimws(get(cl))), .N]
    if (n) err(n, " row(s) with blank '", cl, "'")
}

# 4. recognised platform
bad_pf <- setdiff(unique(na.omit(ss$platform)), PLATFORMS_ALL)
if (length(bad_pf)) err("unrecognised platform(s): ", paste(bad_pf, collapse = ", "),
                        " (allowed: ", paste(PLATFORMS_ALL, collapse = ", "), ")")

ss[, topology := fifelse(platform %in% PLATFORMS_PLATE, "plate", "barcoded")]

# 5. one platform and one group per sample_id
invisible(ss[, {
    if (uniqueN(platform) > 1L) err("sample_id '", .BY$sample_id, "' has multiple platforms: ",
                                    paste(unique(platform), collapse = ", "))
    if (uniqueN(group) > 1L)    err("sample_id '", .BY$sample_id, "' has multiple groups: ",
                                    paste(unique(group), collapse = ", "))
    NULL
}, by = sample_id])

# 3. referenced files. Each distinct path is checked once.
check_file <- function(paths, label, need_index = FALSE, need_2col = FALSE) {
    for (p in unique(na.omit(paths))) {
        if (!file.exists(p)) { err(label, " missing: ", p); next }
        if (file.size(p) == 0) { err(label, " is empty: ", p); next }
        if (need_index && !any(file.exists(paste0(p, ".bai"), sub("\\.bam$", ".bai", p))))
            err("BAM index (.bai) missing for: ", p)
        if (need_2col) {
            ln <- tryCatch(readLines(p, n = 1L, warn = FALSE), error = function(e) character(0))
            if (length(ln) && length(strsplit(ln, "\t", fixed = TRUE)[[1]]) < 2L)
                err(label, " is not 2-column tab-separated: ", p)
        }
    }
}
check_file(ss$bam,             "BAM",             need_index = TRUE)
check_file(ss$celltype_annot,  "celltype_annot",  need_2col  = TRUE)
check_file(ss$called_barcodes, "called_barcodes")

# ── annotation cache ─────────────────────────────────────────────────────────
annot_cache <- new.env(parent = emptyenv())
read_annot <- function(p) {
    key <- digest_key <- p
    if (!is.null(annot_cache[[key]])) return(annot_cache[[key]])
    dt <- tryCatch(
        fread(p, header = FALSE, sep = "\t", colClasses = "character",
              select = 1:2, col.names = c("cell_key", "celltype")),
        error = function(e) NULL
    )
    if (!is.null(dt)) {
        dt[, cell_key := trimws(cell_key)]
        # Cell-type labels become FILENAMES (<group>__<celltype>.bam), so every
        # character outside [A-Za-z0-9._-] collapses to "_". A no-op on all five
        # current datasets (verified), present so an unusual label cannot produce
        # an unquotable path or a name that disagrees with CELLTYPE_SCOPE.tsv.
        dt[, celltype := gsub("[^A-Za-z0-9._-]+", "_", trimws(celltype))]
        dt <- dt[nzchar(cell_key) & nzchar(celltype)]
    }
    annot_cache[[key]] <- dt
    dt
}

# 6. plate-based: cell_id required, unique, and present in its annotation
plate <- ss[topology == "plate"]
if (nrow(plate)) {
    nb <- plate[is.na(cell_id) | !nzchar(trimws(cell_id)), .N]
    if (nb) err(nb, " plate-based row(s) with blank 'cell_id' (required when one BAM = one cell)")
    dup <- plate[!is.na(cell_id), .N, by = cell_id][N > 1L]
    if (nrow(dup)) err(nrow(dup), " duplicated cell_id(s) across the sheet, e.g. ",
                       paste(head(dup$cell_id, 3), collapse = ", "))
    for (ap in unique(plate$celltype_annot)) {
        a <- read_annot(ap); if (is.null(a)) next
        ids <- plate[celltype_annot == ap & !is.na(cell_id), trimws(cell_id)]
        absent <- setdiff(ids, a$cell_key)
        if (length(absent)) err(length(absent), " cell_id(s) absent from ", basename(ap),
                                ", e.g. ", paste(head(absent, 3), collapse = ", "))
    }
}

# 7. barcoded: cell keys unique within a sample_id
bc <- ss[topology == "barcoded"]
if (nrow(bc)) {
    for (sid in unique(bc$sample_id)) {
        aps <- unique(bc[sample_id == sid, celltype_annot])
        keys <- unique(rbindlist(lapply(aps, read_annot), use.names = TRUE, fill = TRUE))
        if (!nrow(keys)) next
        d <- keys[, .N, by = cell_key][N > 1L]
        if (nrow(d)) err("sample_id '", sid, "': ", nrow(d),
                         " cell key(s) appear more than once in its annotation, e.g. ",
                         paste(head(d$cell_key, 3), collapse = ", "))
    }
}

# ── structural report ────────────────────────────────────────────────────────
say("")
say("── structure ─────────────────────────────────────────────")
say("groups      : ", uniqueN(ss$group), "  (", paste(sort(unique(ss$group)), collapse = ", "), ")")
say("samples     : ", uniqueN(ss$sample_id))
say("platform(s) : ", paste(sort(unique(ss$platform)), collapse = ", "),
    "   topology: ", paste(sort(unique(ss$topology)), collapse = ", "))
lanes <- ss[, .N, by = sample_id]
say("lanes/sample: min ", min(lanes$N), "  median ", median(lanes$N), "  max ", max(lanes$N))
if (any(ss$topology == "plate") && any(ss$topology == "barcoded"))
    warn("sheet mixes plate-based and barcoded platforms; the split handles both but they cannot share a group")

say("")
say("── cell types per group (>= ", MIN_CELLS, " cells, pooled across samples) ──")
for (g in sort(unique(ss$group))) {
    sub <- ss[group == g]

    # A cell is identified by (sample_id, cell_key) — which is exactly what the
    # sample_id prefix makes explicit. Keying on cell_key alone would merge cells
    # that different samples happen to share, undercounting the pool.
    per <- rbindlist(lapply(unique(sub$sample_id), function(sid) {
        aps <- unique(sub[sample_id == sid, celltype_annot])
        k <- unique(rbindlist(lapply(aps, read_annot), use.names = TRUE, fill = TRUE))
        if (!nrow(k)) return(NULL)
        if (sub[sample_id == sid, topology[1]] == "plate")
            k <- k[cell_key %in% trimws(na.omit(sub[sample_id == sid, cell_id]))]
        if (nrow(k)) data.table(sample_id = sid, cell_key = k$cell_key, celltype = k$celltype) else NULL
    }), use.names = TRUE, fill = TRUE)
    if (!nrow(per)) { say("  ", g, ": annotation unreadable or empty"); next }

    cts  <- per[, .N, by = celltype][order(-N)]
    keep <- cts[N >= MIN_CELLS]
    say("  ", g, ": ", nrow(keep), " of ", nrow(cts), " cell types kept, ",
        sum(keep$N), " of ", sum(cts$N), " cells")

    # Split plan for this group, taken from the SAME `per`/`keep` objects the
    # scope table is built from — so the BAMs the splitter emits are by
    # construction the ones CELLTYPE_SCOPE.tsv claims were kept.
    if (!is.na(PLAN_DIR)) {
        plan_data[[g]] <- list(
            lanes = sub[, .(sample_id, bam, celltype_annot, platform, topology,
                            cell_id = fifelse(is.na(cell_id), "", trimws(cell_id)))],
            cells = per[celltype %in% keep$celltype][order(sample_id, celltype, cell_key)],
            keep  = sort(keep$celltype)
        )
    }
    if (!nrow(keep)) warn("group '", g, "' has no cell type reaching ", MIN_CELLS,
                          " cells (largest is ", if (nrow(cts)) max(cts$N) else 0, ")")

    # cross-sample key collisions — what the barcode prefix neutralises
    n_coll <- 0L; n_conf <- 0L; pct <- 0
    if (sub$topology[1] == "barcoded" && uniqueN(sub$sample_id) > 1L) {
        coll <- per[, .(n_samples = uniqueN(sample_id), n_types = uniqueN(celltype)),
                    by = cell_key][n_samples > 1L]
        if (nrow(coll)) {
            n_coll <- nrow(coll); n_conf <- nrow(coll[n_types > 1L])
            pct <- round(100 * n_coll / uniqueN(per$cell_key), 1)
            say("      cross-sample key collisions: ", n_coll, " of ", uniqueN(per$cell_key),
                " keys (", pct, "%), ", n_conf, " with conflicting cell types",
                " — neutralised by the sample_id prefix")
            if (pct >= 50)
                warn("group '", g, "': ", pct, "% of cell keys recur across samples — this looks like a ",
                     "fixed plate-barcode set, so (sample_id, barcode) is the only valid cell identity ",
                     "and the prefix is mandatory, not optional")
        }
    }

    scope_rows[[length(scope_rows) + 1L]] <- data.table(
        dataset            = DATASET,
        group              = g,
        samplesheet        = sheet_path,
        platform           = paste(sort(unique(sub$platform)), collapse = "|"),
        topology           = sub$topology[1],
        n_samples          = uniqueN(sub$sample_id),
        n_lanes            = nrow(sub),
        min_cells          = MIN_CELLS,
        n_celltypes_total  = nrow(cts),
        n_celltypes_kept   = nrow(keep),
        n_cells_total      = sum(cts$N),
        n_cells_kept       = sum(keep$N),
        pct_cells_kept     = if (sum(cts$N)) round(100 * sum(keep$N) / sum(cts$N), 2) else NA_real_,
        largest_celltype   = if (nrow(cts)) max(cts$N) else 0L,
        n_keys_colliding   = n_coll,
        pct_keys_colliding = pct,
        n_keys_conflicting = n_conf
    )
}

# ── scope TSV (supplementary table) ──────────────────────────────────────────
# Records the cell-type scope actually applied, per dataset and group, so a
# non-zero --min-cells can never be applied without leaving a trace. Rows for a
# dataset are replaced on re-run so the file stays idempotent.
# Gated on `!length(errs)` for the same reason as the split plan: a sheet that
# does not validate was never applied to anything, so recording its scope in a
# supplementary table would be a false record — and worse, it would REPLACE the
# valid rows of a previous run for that dataset.
if (!is.na(TSV_OUT) && length(scope_rows) && !length(errs)) {
    new <- rbindlist(scope_rows, use.names = TRUE, fill = TRUE)
    new[, validated_utc := format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")]
    if (file.exists(TSV_OUT)) {
        old <- fread(TSV_OUT, colClasses = "character")
        if ("dataset" %in% names(old)) old <- old[dataset != DATASET]
        new <- rbindlist(list(old, new), use.names = TRUE, fill = TRUE)
    }
    dir.create(dirname(TSV_OUT), recursive = TRUE, showWarnings = FALSE)
    fwrite(new, TSV_OUT, sep = "\t")
    say("scope TSV   : ", TSV_OUT, "  (", nrow(new), " rows)")
}

# ── split plan ───────────────────────────────────────────────────────────────
# Written ONLY when the sheet has no errors, so a stale or partial plan can
# never be handed to the splitter. Four files per unit, all tab-separated:
#
#   _expected.txt      sorted unit (group) names, one per line. The array index
#                      IS the line number, which is what lets the submitter size
#                      `-t 1-N` and the driver check completeness against the
#                      same list the tasks read.
#   _plan_<unit>.tsv   one row per LANE: sample_id, bam, celltype_annot,
#                      platform, topology, cell_id.
#   _cells_<unit>.tsv  sample_id, cell_key, celltype — KEPT cell types only,
#                      already sanitised and min-cells filtered.
#   _keep_<unit>.txt   sorted kept cell-type names, one per line.
#
# Everything dataset-specific therefore lives in the samplesheet, and the
# splitter is a pure consumer of these files.
if (!is.na(PLAN_DIR) && !length(errs)) {
    dir.create(PLAN_DIR, recursive = TRUE, showWarnings = FALSE)
    # Stale units from a previous run with different groups would otherwise be
    # picked up by the driver's completeness check.
    unlink(list.files(PLAN_DIR, pattern = "^_(plan|cells|keep)_.*\\.(tsv|txt)$", full.names = TRUE))
    units <- sort(names(plan_data))
    writeLines(units, file.path(PLAN_DIR, "_expected.txt"))
    for (g in units) {
        p <- plan_data[[g]]
        fwrite(p$lanes, file.path(PLAN_DIR, paste0("_plan_",  g, ".tsv")), sep = "\t")
        fwrite(p$cells, file.path(PLAN_DIR, paste0("_cells_", g, ".tsv")), sep = "\t")
        writeLines(p$keep, file.path(PLAN_DIR, paste0("_keep_", g, ".txt")))
    }
    say("split plan  : ", PLAN_DIR, "  (", length(units), " unit(s): ",
        paste(units, collapse = ", "), ")")
}

# ── verdict ──────────────────────────────────────────────────────────────────
say("")
if (length(warns)) {
    cat("WARNINGS (", length(warns), "):\n", sep = "")
    cat(paste0("  - ", warns, collapse = "\n"), "\n", sep = "")
}
if (length(errs)) {
    cat("ERRORS (", length(errs), "):\n", sep = "")
    cat(paste0("  - ", errs, collapse = "\n"), "\n", sep = "")
    if (!is.na(PLAN_DIR)) cat("  (no split plan written)\n")
    cat("\nFAIL: ", sheet_path, "\n", sep = "")
    quit(status = 1L)
}
cat("PASS: ", sheet_path, "\n", sep = "")
quit(status = 0L)
