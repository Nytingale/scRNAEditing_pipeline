#!/usr/bin/env Rscript
# molecule_filter_v8.R <tier>   (Reviewer 2, major 1; user 2026-10-22)
# Main-analysis call sets (published <prefix>_step5_nomm.tsv.gz, NTC tasks of filter_sweep2/manifest.tsv) restricted
# to sites whose genomic variant base is carried by >= MIN_MOL molecules (distinct cell barcode + UMI) in the
# configuration's own BAM (molecule_counts_v8.py), as a sensitivity check of the >= 3 variant reads filter.
# The genomic variant base is the second letter of varType, complemented on the minus strand.
# Input : MW/bams.tsv (tier, unit, bam_type, bam, min_bq, counts), counts tables, published step-5 tables
# Output: MW/results/<rel_dir>/<prefix>_step5_nomm.tsv.gz, MW/summary/<tier>.tsv
# MOL_DIR (default gt_eval_v8/molecules) and MOL_MANIFEST (default filter_sweep2/manifest.tsv; an optional calls
# column gives the step-5 table, else PM/<rel_dir>/<prefix>_step5_nomm.tsv.gz) for the v9 run (molecule_check_v9.sh).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
TIER <- args[1]
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
MW <- Sys.getenv("MOL_DIR", file.path(PM, "gt_eval_v8", "molecules"))
MIN_MOL <- as.integer(Sys.getenv("MIN_MOL", "3"))
COMP <- c(A = "T", C = "G", G = "C", T = "A")

man <- fread(Sys.getenv("MOL_MANIFEST", file.path(PM, "gt_eval_v8", "filter_sweep2", "manifest.tsv")))[tier == TIER]
man[, c("unit", "bam_type") := tstrsplit(rel_dir, "/", keep = c(2L, 5L))]
if (!"calls" %in% names(man)) man[, calls := file.path(PM, rel_dir, paste0(prefix, "_step5_nomm.tsv.gz"))]
man <- man[bam_type != "sample"]
bams <- fread(file.path(MW, "bams.tsv"))[tier == TIER]
stopifnot(nrow(man) > 0, all(paste(man$unit, man$bam_type) %in% paste(bams$unit, bams$bam_type)))

cnt_cache <- list()
counts_for <- function(u, bt) {
  k <- paste(u, bt)
  if (is.null(cnt_cache[[k]])) {
    f <- bams[unit == u & bam_type == bt, counts]
    cnt_cache[[k]] <<- fread(f, select = c("chr", "pos", "ref", paste0("n_", names(COMP)), paste0("m_", names(COMP))))
  }
  cnt_cache[[k]]
}

summ <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  r <- man[i]
  f <- r$calls
  d <- fread(f)
  out_dir <- file.path(MW, "results", r$rel_dir); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, paste0(r$prefix, "_step5_nomm.tsv.gz"))
  if (nrow(d) == 0) { fwrite(d, out, sep = "\t"); return(data.table(id = r$id, n = 0L)) }
  stopifnot(all(d$strand %in% c("+", "-")), all(nchar(d$varType) == 2))
  g_ref <- substr(d$varType, 1, 1); g_alt <- substr(d$varType, 2, 2)
  m <- d$strand == "-"
  g_ref[m] <- COMP[g_ref[m]]; g_alt[m] <- COMP[g_alt[m]]
  cn <- counts_for(r$unit, r$bam_type)
  x <- cn[data.table(chr = d$chr, pos = d$pos), on = c("chr", "pos")]
  stopifnot(!anyNA(x$ref))
  idx <- cbind(seq_len(nrow(d)), match(g_alt, names(COMP)))
  alt_reads <- as.matrix(x[, paste0("n_", names(COMP)), with = FALSE])[idx]
  alt_mol <- as.matrix(x[, paste0("m_", names(COMP)), with = FALSE])[idx]
  keep <- alt_mol >= MIN_MOL
  fwrite(d[keep], out, sep = "\t")
  ag <- d$varType == "AG"
  data.table(id = r$id, caller = r$caller, unit = r$unit, bam_type = r$bam_type, prefix = r$prefix, n = nrow(d),
             n_kept = sum(keep), n_ag = sum(ag), n_ag_kept = sum(keep & ag), ref_mismatch = sum(x$ref != g_ref),
             alt_reads = sum(alt_reads), alt_mol = sum(alt_mol), alt_reads_ag = sum(alt_reads[ag]),
             alt_mol_ag = sum(alt_mol[ag]), n_alt_reads_lt3 = sum(alt_reads < 3))
}), fill = TRUE)
dir.create(file.path(MW, "summary"), showWarnings = FALSE)
fwrite(summ, file.path(MW, "summary", paste0(TIER, ".tsv")), sep = "\t")
cat(TIER, ":", nrow(summ), "call sets; kept", round(100 * sum(summ$n_kept) / sum(summ$n), 1), "% of sites,",
    round(100 * sum(summ$n_ag_kept, na.rm = TRUE) / sum(summ$n_ag, na.rm = TRUE), 1), "% of A-to-G; ref mismatches",
    sum(summ$ref_mismatch, na.rm = TRUE), "\n")
