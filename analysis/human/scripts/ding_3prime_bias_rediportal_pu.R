#!/usr/bin/env Rscript
# =============================================================================
# ding_3prime_bias_rediportal_pu.R   (2026-10-21; Reviewer 3, major 2: editing sites missed because of 3' bias)
# ding_3prime_bias_rediportal.R on the genome-wide pileups of ding_3prime_bias_pileup.sh, so that both platforms come
# from the same donor (Smart-seq2 PBMC1 only, as the 10x data) and each can be read after its own preprocessing:
# 10x 3' v3 UMI-deduplicated (x10u) vs Smart-seq2 PBMC1 all reads (ss2p1) or per-cell coordinate-deduplicated
# (ss2p1d). Without UMIs, coordinate deduplication of full-length Smart-seq2 collapses distinct molecules sharing a
# start in highly expressed genes, so all reads (ss2p1) is the main Smart-seq2 input and ss2p1d a sensitivity check.
# Depth (DEPTH):
#   sequenced  equal sequencing effort: the platform with more sequenced (primary, before deduplication) reads is
#              binomially thinned by the ratio of sequenced reads (10x 473.8 M vs Smart-seq2 PBMC1 426.7 M)
#   processed  equal reads after each platform's preprocessing (read_counts.tsv of the BAM each pileup came from)
#   subsampled X10=x10us20261024 SS2=ss2p1ds<seed>: the pipeline's deduplicated BAMs subsampled to exactly
#              215,322,584 reads each (Smart-seq2 as whole proper pairs; ding_ss2_dedup_subsample.sh), so no thinning
#   covered : coverage >= 10;   edited : coverage >= 10 and >= 3 edited reads (as in the other analyses) (REDIportal Ed base on the site's strand;
#   pileups are unstranded and keep A/T reference positions with >= 10 reads)
# Region (REDIportal Func.refGene) and genomic distance to the nearest RefSeq transcript 3' end on the site's strand.
# Thinning: each site's edited and other read counts are binomially thinned (seed SEED, written to the outputs).
# Usage : X10=x10u SS2=ss2p1 DEPTH=sequenced Rscript ding_3prime_bias_rediportal_pu.R
# Output: human/jzlxv3/three_prime_bias/rp_<X10>_<SS2>_<DEPTH>_{sites.tsv.gz, region.tsv, tes.tsv}
# =============================================================================
suppressPackageStartupMessages(library(data.table))
H <- file.path(Sys.getenv("SCRNA_DIR", "."), "human")
OUT <- file.path(H, "jzlxv3", "three_prime_bias")
PU <- file.path(OUT, "pileup")
MIN_COV <- 10L
MIN_ALT <- 3L
SEED <- 20261021L
set.seed(SEED)
X10 <- Sys.getenv("X10", "x10u"); SS2 <- Sys.getenv("SS2", "ss2p1"); DEPTH <- Sys.getenv("DEPTH", "sequenced")
SUB_SEEDS <- c(20261024L, 20261025L, 20261026L)
stopifnot(X10 %in% c("x10", "x10u", "x10us20261024"), SS2 %in% c("ss2p1", "ss2p1d", paste0("ss2p1ds", SUB_SEEDS)),
          DEPTH %in% c("sequenced", "processed", "subsampled"),
          (DEPTH == "subsampled") == grepl("^ss2p1ds", SS2), (DEPTH == "subsampled") == grepl("^x10us", X10))
PREFIX <- paste("rp", X10, SS2, DEPTH, sep = "_")

PLAT_BAM <- c(x10 = "jzlxv3/wasp_refbias/merged_all.bam", x10u = "jzlxv3/results_v8/PBMC1/02_dedup/PBMC1_umiDup.bam",
              ss2p1 = "jzlss2/wasp_refbias/merged_pbmc1.bam",
              ss2p1d = "jzlss2/results_v8/PBMC1/02_dedup/PBMC1_rmDup.bam")
PLAT_BAM[paste0("ss2p1ds", SUB_SEEDS)] <- sprintf("jzlss2/results_v8/PBMC1/02_dedup/subsample/PBMC1_rmDup_sub%d.bam",
                                                  SUB_SEEDS)
SEQ_BAM <- c(x10 = "x10", x10u = "x10", ss2p1 = "ss2p1", ss2p1d = "ss2p1")
rc <- fread(file.path(OUT, "read_counts.tsv"))
n_reads <- function(p) rc[bam == file.path(H, PLAT_BAM[[p]]), reads_primary]
if (DEPTH == "subsampled") {
  n_10x <- fread(file.path(H, "jzlxv3/results_v8/PBMC1/02_dedup/subsample",
                           sprintf("PBMC1_umiDup_sub%s.count", sub("x10us", "", X10))))$reads_primary
  n_ss2 <- fread(file.path(H, "jzlss2/results_v8/PBMC1/02_dedup/subsample",
                           sprintf("PBMC1_rmDup_sub%s.count", sub("ss2p1ds", "", SS2))))$reads_primary
  stopifnot(n_10x == n_ss2)
} else {
  n_10x <- n_reads(if (DEPTH == "sequenced") SEQ_BAM[[X10]] else X10)
  n_ss2 <- n_reads(if (DEPTH == "sequenced") SEQ_BAM[[SS2]] else SS2)
}
stopifnot(length(n_10x) == 1, length(n_ss2) == 1)
THIN <- if (n_ss2 > n_10x) "ss2" else if (n_10x > n_ss2) "x10" else "none"
f_thin <- min(n_10x, n_ss2) / max(n_10x, n_ss2)
message(sprintf("%s: 10x %s %.1f M, Smart-seq2 %s %.1f M (%s reads): thinning %s by %.3f", PREFIX, X10, n_10x / 1e6,
                SS2, n_ss2 / 1e6, DEPTH, THIN, f_thin))

# TABLE1 has 25 header names but 22 data columns (the three GeneDetail columns are absent): Func.refGene is column 14
rp <- fread(cmd = paste("zcat", file.path(H, "ref_annots", "TABLE1_hg38_v2.txt.gz"), "| cut -f1-5,14"))
setnames(rp, c("chr", "pos", "ref", "ed", "strand", "func"))
CHRS <- intersect(paste0("chr", c(1:22, "X")), unique(rp$chr))

read_pu <- function(plat, cc, thin) {
  p <- fread(file.path(PU, sprintf("%s_%s.tsv.gz", plat, cc)), colClasses = list(character = "ref"))
  for (b in c("A", "C", "G", "T")) set(p, j = b, value = p[[b]] + p[[tolower(b)]])
  m <- merge(rp[chr == cc, .(chr, pos, ed)], p[, .(chr, pos, coverage, A, C, G, T)], by = c("chr", "pos"))
  m[, alt := fcase(ed == "G", G, ed == "C", C, ed == "A", A, ed == "T", T)]
  if (thin) {
    rest <- m$coverage - m$alt
    m[, alt := rbinom(.N, alt, f_thin)]
    m[, coverage := alt + rbinom(.N, rest, f_thin)]
  }
  m[, .(chr, pos, ed, cov = coverage, alt)]
}
pu <- rbindlist(lapply(CHRS, function(cc) {
  x <- read_pu(X10, cc, THIN == "x10")
  s <- read_pu(SS2, cc, THIN == "ss2")
  merge(x[, .(chr, pos, ed, cov_10x = cov, alt_10x = alt)], s[, .(chr, pos, ed, cov_ss2 = cov, alt_ss2 = alt)],
        by = c("chr", "pos", "ed"), all = TRUE)
}))
d <- merge(rp[chr %chin% CHRS], pu, by = c("chr", "pos", "ed"), all.x = TRUE)
for (v in c("cov_10x", "alt_10x", "cov_ss2", "alt_ss2")) set(d, which(is.na(d[[v]])), v, 0L)
d[, region := fcase(func == "UTR3", "3' UTR", func == "exonic", "Coding exon", func == "UTR5", "5' UTR",
                    func == "intronic", "Intron", grepl("^ncRNA", func), "Noncoding RNA",
                    func %chin% c("downstream", "upstream", "upstream;downstream"), "Up- or downstream",
                    default = "Intergenic and other")]

gtf <- fread(cmd = paste("awk -F'\\t' '$3==\"transcript\"'", file.path(H, "ref_annots", "hg38.ncbiRefSeq.gtf")),
             header = FALSE, select = c(1, 4, 5, 7))
setnames(gtf, c("chr", "start", "end", "strand"))
tes_dt <- unique(gtf[, .(chr, strand, tes = fifelse(strand == "+", end, start))])[, tes_pos := tes]
setkey(tes_dt, chr, strand, tes)
d[, tes := pos]
d[, d_tes := abs(tes_dt[d, on = .(chr, strand, tes), roll = "nearest", tes_pos] - pos)]
d[, tes := NULL]
d[, `:=`(covered_10x = cov_10x >= MIN_COV, covered_ss2 = cov_ss2 >= MIN_COV,
         edited_10x = cov_10x >= MIN_COV & alt_10x >= MIN_ALT, edited_ss2 = cov_ss2 >= MIN_COV & alt_ss2 >= MIN_ALT)]
fwrite(d[covered_10x | covered_ss2], file.path(OUT, paste0(PREFIX, "_sites.tsv.gz")), sep = "\t")

cls <- function(a, b) fcase(a & b, "Both", a, "10x 3' v3 only", b, "Smart-seq2 only", default = "Neither")
reg <- rbind(
  d[, .(what = "covered", class = cls(covered_10x, covered_ss2), region)],
  d[, .(what = "edited", class = cls(edited_10x, edited_ss2), region)])[class != "Neither", .N, by = .(what, region, class)]
reg <- rbind(reg, reg[, .(region = "All", N = sum(N)), by = .(what, class)])
reg[, pct := 100 * N / sum(N), by = .(what, region)]
reg[, `:=`(x10 = X10, ss2 = SS2, depth = DEPTH, f_thin = f_thin, thinned = THIN, seed = SEED)]
fwrite(reg, file.path(OUT, paste0(PREFIX, "_region.tsv")), sep = "\t")
print(dcast(reg, what + region ~ class, value.var = "pct", fun.aggregate = function(z) round(z, 1)))

BRK <- c(-1, 500, 1000, 2000, 5000, 10000, 20000, 50000, Inf)
LAB <- c("0-0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", "20-50", ">50")
d[, tes_bin := factor(LAB[findInterval(d_tes, BRK, left.open = TRUE)], LAB)]
EXONIC <- c("3' UTR", "Coding exon", "5' UTR")
tb <- rbind(d[!is.na(tes_bin)][, sites := "All"], d[!is.na(tes_bin) & region %chin% EXONIC][, sites := "Exonic"])[
  , .(n_sites = .N, covered_10x = sum(covered_10x), covered_ss2 = sum(covered_ss2), edited_10x = sum(edited_10x),
      edited_ss2 = sum(edited_ss2)), by = .(sites, tes_bin)][order(sites, tes_bin)]
tb[, `:=`(x10 = X10, ss2 = SS2, depth = DEPTH, f_thin = f_thin, thinned = THIN, seed = SEED)]
fwrite(tb, file.path(OUT, paste0(PREFIX, "_tes.tsv")), sep = "\t")
print(tb)
message(sprintf("%s: %d REDIportal sites on %d chromosomes", PREFIX, nrow(d), length(CHRS)))
