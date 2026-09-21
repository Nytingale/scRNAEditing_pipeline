#!/usr/bin/env Rscript
# ============================================================================
# test_annot_chunk_equivalence.R
# Proves the chunked STEP 5 annotation lookup in downstream_filter.R returns
# EXACTLY the same voted table as the single-block version it replaced.
# ============================================================================
# The chunking exists because the single-block lookup OOM-killed at 128 GB on
# 31.2 M sites (Replogle merged NTC / reditools3), and because findOverlaps
# expands to one row per (site, feature) pair, peak memory scaled with site
# count. The equivalence argument is that processAnnotatedBED() aggregates
# strictly by (chr, pos) and every chunk is matched against the FULL annotation,
# so a site's feature set is never split. This script tests that claim on real
# annotation data rather than asserting it.
#
# Usage: Rscript test_annot_chunk_equivalence.R <strand_annots.RData> [N_SITES]
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
ANNOT <- if (length(args) >= 1) args[1] else
  "/illumina-isi03/scratch/deep_learning/yng7/scRNAEditing/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData"
N_SITES <- if (length(args) >= 2) as.integer(args[2]) else 300000L

stopifnot(file.exists(ANNOT))
cat("annotation :", ANNOT, "\n")
load(ANNOT)   # plus.annot, minus.annot
cat("plus.annot :", length(plus.annot), " features\n")
cat("minus.annot:", length(minus.annot), " features\n")

# ── the production voting function, copied verbatim ─────────────────────────
processAnnotatedBED = function(df) {
  colnames(df) = c("chr", "pos", "gene_name", "strand", "gene_type")
  df$gene_type[is.na(df$gene_type)] = "NA"
  type_names = c("NA", "gene", "transcript", "intergenic", "intron",
                 "5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR")
  ranks = c(0, 0, 0, 1, 2, rep(3, 6))
  map = setNames(ranks, type_names)
  df$ranktype = map[df$gene_type]
  df$ranktype[is.na(df$ranktype)] = 0
  df = as.data.table(df)
  df = df[, .(genetype  = paste(unique(gene_type), collapse = ","),
              genename  = paste(unique(gene_name), collapse = ","),
              ranktype  = max(ranktype)),
          by = .(chr, pos)]
  df$genetype = stringr::str_replace_na(df$genetype)
  df$genename = stringr::str_replace_na(df$genename)
  as.data.frame(df)
}

vote_overlaps = function(gr, annot) {
  hits = findOverlaps(gr, annot)
  if (length(hits) == 0L) {
    return(data.frame(chr = character(0), pos = numeric(0), genetype = character(0),
                      genename = character(0), ranktype = numeric(0),
                      stringsAsFactors = FALSE))
  }
  matched = gr[queryHits(hits)]
  mcols(matched) = cbind.data.frame(mcols(matched), mcols(annot[subjectHits(hits)]))
  d = as.data.frame(matched)[, c("seqnames", "end", "id", "strand", "geneType")]
  colnames(d) = c("chr", "pos", "gene_name", "strand", "gene_type")
  d$chr = as.character(d$chr)
  processAnnotatedBED(unique(d))
}

# ── a realistic site set: real annotated positions plus unannotated gaps ───
set.seed(42)
src <- c(plus.annot, minus.annot)
idx <- sample(seq_along(src), min(N_SITES, length(src)))
site_chr <- as.character(seqnames(src)[idx])
site_pos <- start(src)[idx] + sample(0:400, length(idx), replace = TRUE)
annot_df <- unique(data.frame(chr = site_chr, pos = site_pos, stringsAsFactors = FALSE))
annot_df$site <- paste0(annot_df$chr, ":", annot_df$pos)
annot_df <- annot_df[!duplicated(annot_df$site), ]
cat("test sites :", nrow(annot_df), "\n")

run_blocks <- function(chunk) {
  n <- nrow(annot_df)
  starts <- seq(1L, n, by = chunk)
  pp <- vector("list", length(starts)); mm <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    rows <- starts[k]:min(starts[k] + chunk - 1L, n)
    pk <- as.numeric(annot_df$pos[rows])
    g <- GRanges(annot_df$chr[rows], IRanges(pk, pk))
    pp[[k]] <- vote_overlaps(g, plus.annot)
    mm[[k]] <- vote_overlaps(g, minus.annot)
  }
  list(plus  = as.data.frame(rbindlist(pp, use.names = TRUE, fill = TRUE)),
       minus = as.data.frame(rbindlist(mm, use.names = TRUE, fill = TRUE)),
       n_blocks = length(starts))
}

norm <- function(d) {
  if (!nrow(d)) return(d)
  d$chr <- as.character(d$chr); d$pos <- as.numeric(d$pos)
  d <- d[order(d$chr, d$pos), c("chr", "pos", "genetype", "genename", "ranktype")]
  rownames(d) <- NULL
  d
}

cat("\n-- single block (the ORIGINAL code path) --\n")
one <- run_blocks(nrow(annot_df))
cat("   blocks:", one$n_blocks, "  plus rows:", nrow(one$plus),
    "  minus rows:", nrow(one$minus), "\n")

ok <- TRUE
for (chunk in c(100000L, 25000L, 7000L)) {
  many <- run_blocks(chunk)
  p_eq <- identical(norm(one$plus),  norm(many$plus))
  m_eq <- identical(norm(one$minus), norm(many$minus))
  cat(sprintf("-- chunk %6d -> %3d blocks : plus %s, minus %s\n",
              chunk, many$n_blocks,
              if (p_eq) "IDENTICAL" else "DIFFERS",
              if (m_eq) "IDENTICAL" else "DIFFERS"))
  if (!p_eq || !m_eq) {
    ok <- FALSE
    a <- norm(one$plus); b <- norm(many$plus)
    cat("   plus  rows:", nrow(a), "vs", nrow(b), "\n")
    a <- norm(one$minus); b <- norm(many$minus)
    cat("   minus rows:", nrow(a), "vs", nrow(b), "\n")
  }
}

cat("\n")
if (ok) {
  cat("PASS: chunked voting is identical to single-block voting at every chunk size\n")
  quit(status = 0)
} else {
  cat("FAIL: chunking changed the result\n")
  quit(status = 1)
}
