#!/usr/bin/env Rscript
# =============================================================================
# ding_3prime_bias_denovo.R   (2026-10-18; Reviewer 3, major 2: editing sites missed because of 3' bias)
# De novo (REDIportal-independent) version of ding_3prime_bias_rediportal.R. Candidate sites are called in the same
# way on both platforms from genome-wide pileups of Ding et al. human PBMC 10x 3' v3 (PBMC1) and Smart-seq2
# (default the 384 PBMC1 cells, same donor; SS2_PLATFORM=ss2 for PBMC1 + PBMC2, whose PBMC2-private SNPs at ~25%
# dominate the Smart-seq2 candidates). The platform with more primary reads is binomially thinned to the other's
# (read_counts.tsv from ding_ss2_pbmc1_merge.sh; PBMC1 only: 10x 473.8 M thinned to 426.7 M).
#   candidate : coverage >= 10, >= MIN_ALT reads of the mismatch (default 3, as in the other analyses; the 10-18 / 10-19
#               runs used 2), mismatch fraction < 0.9, not a dbSNP156 common position
# Both pileups are unstranded and keep A/T reference positions only, so each mismatch type is read on the strand on
# which the reference base is A: A-to-G (ref A > G, + strand; ref T > C, - strand), and the controls A-to-C (A > C /
# T > G) and A-to-T (A > T / T > A). The A-to-G share of all candidates is a caller-independent quality measure.
# Region (RefSeq, on that strand; CDS > 3' UTR > 5' UTR > other exon > intron > intergenic) and genomic distance to the
# nearest RefSeq transcript 3' end on that strand.
# Inputs : human/jzlxv3/three_prime_bias/pileup/{x10,x10u,ss2p1,ss2,ss2p1d}_chrN.tsv.gz (ding_3prime_bias_pileup.sh),
#          human/jzlxv3/three_prime_bias/read_counts.tsv (ding_ss2_pbmc1_merge.sh; dedup BAMs ding_dedup_read_counts.sh),
#          human/pm/gt_eval_v6/snp/dbSNP156_common.bed, ref_annots/{hg38.ncbiRefSeq.gtf, TABLE1_hg38_v2.txt.gz}
# Output : human/jzlxv3/three_prime_bias/denovo_{sites.tsv.gz, region.tsv, tes.tsv, mismatch.tsv}
# =============================================================================
suppressPackageStartupMessages(library(data.table))
H <- file.path(Sys.getenv("SCRNA_DIR", "."), "human")
OUT <- file.path(H, "jzlxv3", "three_prime_bias")
PU <- file.path(OUT, "pileup")
MIN_COV <- 10L
MIN_ALT <- as.integer(Sys.getenv("MIN_ALT", "3"))
MAX_FRAC <- 0.9
CHRS <- paste0("chr", c(1:22, "X"))
PFX <- "denovo_"
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  CHRS <- args
  PFX <- "denovo_test_"
}
SEED <- 20261018L
set.seed(SEED)

# Smart-seq2 pileup: ss2p1 = PBMC1 cells only (same donor as the 10x data; default), ss2 = PBMC1 + PBMC2,
# ss2p1d = PBMC1 per-cell coordinate-deduplicated. 10x pileup: x10 = all reads (default), x10u = UMI-deduplicated.
# Deduplicated pair: X10_PLATFORM=x10u SS2_PLATFORM=ss2p1d (outputs prefixed denovo_dedup_)
SS2 <- Sys.getenv("SS2_PLATFORM", "ss2p1")
X10 <- Sys.getenv("X10_PLATFORM", "x10")
PLAT_BAM <- c(x10 = "jzlxv3/wasp_refbias/merged_all.bam", x10u = "jzlxv3/results_v8/PBMC1/02_dedup/PBMC1_umiDup.bam",
              ss2p1 = "jzlss2/wasp_refbias/merged_pbmc1.bam", ss2 = "jzlss2/wasp_refbias/merged_all.bam",
              ss2p1d = "jzlss2/results_v8/PBMC1/02_dedup/PBMC1_rmDup.bam")
stopifnot(X10 %in% c("x10", "x10u"), SS2 %in% c("ss2p1", "ss2", "ss2p1d"))
# DEPTH: processed = reads of the BAM each pileup came from (default; the 10-18 / 10-19 runs); sequenced = reads
# before deduplication (x10u -> x10, ss2p1d -> ss2p1), i.e. equal sequencing effort
DEPTH <- Sys.getenv("DEPTH", "processed")
stopifnot(DEPTH %in% c("processed", "sequenced"))
SEQ_BAM <- c(x10 = "x10", x10u = "x10", ss2p1 = "ss2p1", ss2 = "ss2", ss2p1d = "ss2p1")
if (MIN_ALT != 2L || DEPTH != "processed") {
  PFX <- sprintf("%s%s_%s_%s_alt%d_", PFX, X10, SS2, DEPTH, MIN_ALT)
} else if (X10 == "x10u" && SS2 == "ss2p1d") {
  PFX <- paste0(PFX, "dedup_")
} else {
  if (X10 != "x10") PFX <- paste0(PFX, X10, "_")
  if (SS2 != "ss2p1") PFX <- paste0(PFX, SS2, "_")
}
rc <- fread(file.path(OUT, "read_counts.tsv"))
n_10x <- rc[endsWith(bam, PLAT_BAM[[if (DEPTH == "sequenced") SEQ_BAM[[X10]] else X10]]), reads_primary]
n_ss2 <- rc[endsWith(bam, PLAT_BAM[[if (DEPTH == "sequenced") SEQ_BAM[[SS2]] else SS2]]), reads_primary]
stopifnot(length(n_10x) == 1, length(n_ss2) == 1)
# the platform with more primary reads is binomially thinned to the other's read number
THIN <- if (n_ss2 > n_10x) "ss2" else "x10"
f_thin <- min(n_10x, n_ss2) / max(n_10x, n_ss2)
message(sprintf("Smart-seq2 = %s (%d reads), 10x = %d reads: thinning %s by %.3f", SS2, n_ss2, n_10x, THIN, f_thin))

# mismatch types on the A strand: ref A -> (G, C, T); ref T -> complements (C, G, A)
MM <- data.table(mm = c("A-to-G", "A-to-C", "A-to-T"), alt_a = c("G", "C", "T"), alt_t = c("C", "G", "A"))

read_pu <- function(plat, chr, thin) {
  f <- file.path(PU, sprintf("%s_%s.tsv.gz", plat, chr))
  stopifnot(file.exists(f))
  p <- fread(f, colClasses = list(character = "ref"))
  for (b in c("A", "C", "G", "T")) set(p, j = b, value = p[[b]] + p[[tolower(b)]])
  p[, rep := ref %chin% c("a", "t")]
  p[, ref := toupper(ref)]
  if (thin) {
    p[, other := pmax(coverage - (A + C + G + T), 0L)]
    for (b in c("A", "C", "G", "T", "other")) set(p, j = b, value = rbinom(nrow(p), p[[b]], f_thin))
    p[, coverage := other + A + C + G + T]
  }
  r <- rbindlist(lapply(seq_len(nrow(MM)), function(k) {
    p[, .(chr, pos, ref, rep, mm = MM$mm[k], cov = coverage,
          alt = fifelse(ref == "A", get(MM$alt_a[k]), get(MM$alt_t[k])))]
  }))
  r[, cand := cov >= MIN_COV & alt >= MIN_ALT & alt / cov < MAX_FRAC]
  r
}

snp <- fread(file.path(H, "pm", "gt_eval_v6", "snp", "dbSNP156_common.bed"), header = FALSE, select = c(1, 3),
             col.names = c("chr", "pos"))
setkey(snp, chr, pos)

d <- rbindlist(lapply(CHRS, function(cc) {
  x <- read_pu(X10, cc, THIN == "x10")
  s <- read_pu(SS2, cc, THIN == "ss2")
  m <- merge(x[, .(chr, pos, ref, rep, mm, cov_10x = cov, alt_10x = alt, cand_10x = cand)],
             s[, .(chr, pos, ref, rep, mm, cov_ss2 = cov, alt_ss2 = alt, cand_ss2 = cand)],
             by = c("chr", "pos", "ref", "rep", "mm"), all = TRUE)
  m <- m[cand_10x %in% TRUE | cand_ss2 %in% TRUE]
  for (v in c("cov_10x", "alt_10x", "cov_ss2", "alt_ss2")) set(m, which(is.na(m[[v]])), v, 0L)
  for (v in c("cand_10x", "cand_ss2")) set(m, which(is.na(m[[v]])), v, FALSE)
  m <- m[!snp[chr == cc], on = .(chr, pos)]
  message(sprintf("%s: %d candidate site x mismatch rows", cc, nrow(m)))
  m
}))
d[, strand := fifelse(ref == "A", "+", "-")]

rp <- fread(cmd = paste("zcat", file.path(H, "ref_annots", "TABLE1_hg38_v2.txt.gz"), "| cut -f1,2"))
setnames(rp, c("chr", "pos"))
d[, rediportal := FALSE][rp, rediportal := TRUE, on = .(chr, pos)]

gtf <- fread(cmd = paste("awk -F'\\t' '$3==\"transcript\" || $3==\"exon\" || $3==\"CDS\" || $3==\"3UTR\" || $3==\"5UTR\"'",
                         file.path(H, "ref_annots", "hg38.ncbiRefSeq.gtf")),
             header = FALSE, select = c(1, 3, 4, 5, 7), col.names = c("chr", "feature", "start", "end", "strand"))
PRIO <- c(CDS = 1L, `3UTR` = 2L, `5UTR` = 3L, exon = 4L, transcript = 5L)
REG <- c("Coding exon", "3' UTR", "5' UTR", "Noncoding exon", "Intron")
gtf[, prio := PRIO[feature]]
gtf <- unique(gtf[, .(chr, strand, start, end, prio)])
setkey(gtf, chr, strand, start, end)
q <- unique(d[, .(chr, strand, start = pos, end = pos)])
ov <- foverlaps(q, gtf, by.x = c("chr", "strand", "start", "end"), nomatch = NULL)[
  , .(prio = min(prio)), by = .(chr, strand, pos = i.start)]
d[ov, prio := i.prio, on = .(chr, strand, pos)]
d[, region := fifelse(is.na(prio), "Intergenic", REG[prio])][, prio := NULL]

tes_dt <- unique(gtf[prio == 5L, .(chr, strand, tes = fifelse(strand == "+", end, start))])[, tes_pos := tes]
setkey(tes_dt, chr, strand, tes)
d[, tes := pos]
d[, d_tes := abs(tes_dt[d, on = .(chr, strand, tes), roll = "nearest", tes_pos] - pos)]
d[, tes := NULL]
fwrite(d, file.path(OUT, paste0(PFX, "sites.tsv.gz")), sep = "\t")

cls <- function(a, b) fcase(a & b, "Both", a, "10x 3' v3 only", b, "Smart-seq2 only")
d[, class := cls(cand_10x, cand_ss2)]

mmt <- rbind(d[, .N, by = .(mm, class)], d[, .(class = "All", N = .N), by = mm])
mmt <- rbind(mmt, d[cand_10x == TRUE, .(class = "10x 3' v3 (all)", .N), by = mm],
             d[cand_ss2 == TRUE, .(class = "Smart-seq2 (all)", .N), by = mm])
mmt[, pct := 100 * N / sum(N), by = class]
fwrite(mmt, file.path(OUT, paste0(PFX, "mismatch.tsv")), sep = "\t")
print(dcast(mmt, class ~ mm, value.var = "pct"))

ag <- d[mm == "A-to-G"]
reg <- ag[, .(N = .N, rediportal = sum(rediportal), alu_soft_masked = sum(rep)), by = .(region, class)]
reg <- rbind(reg, reg[, .(region = "All", N = sum(N), rediportal = sum(rediportal), alu_soft_masked = sum(alu_soft_masked)),
                      by = class])
ctl <- d[mm != "A-to-G", .(n_control = .N), by = .(region, class)]
ctl <- rbind(ctl, ctl[, .(region = "All", n_control = sum(n_control)), by = class])
reg <- merge(reg, ctl, by = c("region", "class"), all.x = TRUE)
reg[is.na(n_control), n_control := 0L]
reg[, `:=`(pct = 100 * N / sum(N), ag_share = 100 * N / (N + n_control)), by = region]
setnames(reg, "alu_soft_masked", "repeat_soft_masked")
fwrite(reg, file.path(OUT, paste0(PFX, "region.tsv")), sep = "\t")
print(reg[order(region, class)])

BRK <- c(-1, 500, 1000, 2000, 5000, 10000, 20000, 50000, Inf)
LAB <- c("0-0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", "20-50", ">50")
d[, tes_bin := factor(LAB[findInterval(d_tes, BRK, left.open = TRUE)], LAB)]
EXONIC <- c("Coding exon", "3' UTR", "5' UTR", "Noncoding exon")
tb <- rbind(d[!is.na(tes_bin)][, sites := "All"], d[!is.na(tes_bin) & region %chin% EXONIC][, sites := "Exonic"])[
  , .(cand_10x = sum(cand_10x), cand_ss2 = sum(cand_ss2), both = sum(cand_10x & cand_ss2)),
  by = .(sites, mm, tes_bin)][order(sites, mm, tes_bin)]
fwrite(tb, file.path(OUT, paste0(PFX, "tes.tsv")), sep = "\t")
print(tb[mm == "A-to-G"])
message(sprintf("seed %d, MIN_ALT %d, DEPTH %s; %s thinned by %.3f; %d candidate rows (%d A-to-G)", SEED, MIN_ALT, DEPTH, THIN, f_thin, nrow(d), nrow(ag)))
