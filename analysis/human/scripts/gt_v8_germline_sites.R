#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_germline_sites.R   (2026-10-13)
# Germline negatives for the germline false-positive rate (user 2026-10-13): PGP1 germline SNVs (donor of
# the iPSC line; whole-genome calls independent of the RNA data) that look like A-to-I editing in the
# NTC reads. A site is a (position, read orientation) with
#   PGP1 SNV         FILTER PASS, genomic REF>ALT = A>G (seen on + reads) or T>C (seen on - reads),
#                    i.e. RNA-frame A>G in 10x read orientation, the same frame as the truth sites
#   not dbSNP common (spectrum12 positions_noSNP: candidate positions with dbSNP156 common removed)
#   NTC evidence     merged NTC (all NTC cells, UMI dedup, BQ/MQ 30) on that orientation: coverage >= 10,
#                    alt >= 3, alt / coverage >= 5% (truth-like: frozen min_edit_reads / min_ctrl_edit)
#   not multiallelic third allele on that orientation < 3 reads or < 5% of depth (truth rule)
# v8 truth positions are disjoint by construction (PGP1 positions are removed from v8).
# Input : pgp1_snps/PGP1_hu43860C_hg38_snvs.vcf.gz, gt_eval_v6/spectrum12/{positions_noSNP.bed,pileups/}
# Output: gt_eval_v8/germline/
#   germline_sites.tsv.gz   set, site, chr, pos, genomic_ref, genomic_alt, read_strand, gt, ntc_cov, ntc_alt,
#                           ntc_vaf, ko1_vaf, ko2_vaf (merged KO arms, same orientation; NA if cov < 10)
#   germline_counts.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
SP <- file.path(PM, "gt_eval_v6", "spectrum12")
OUT <- file.path(PM, "gt_eval_v8", "germline"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fr <- fread(file.path(PM, "gt_eval_v6", "FROZEN_RULES.tsv")); P <- setNames(fr$value, fr$parameter)
K <- as.numeric(P[["min_edit_reads"]]); E <- 0.05; MIN_COV <- as.numeric(P[["min_cov"]])
MIN_THIRD_READS <- 3; MIN_THIRD_FRAC <- 0.05
message("rule: cov >= ", MIN_COV, ", alt >= ", K, ", VAF >= ", E)

v <- fread(cmd = paste("zcat", file.path(PM, "pgp1_snps", "PGP1_hu43860C_hg38_snvs.vcf.gz"), "| grep -v '^##'"),
           select = c(1, 2, 4, 5, 7, 10), col.names = c("chr", "pos", "ref", "alt", "filter", "gt"))
n_vcf <- nrow(v)
v <- v[filter == "PASS" & ((ref == "A" & alt == "G") | (ref == "T" & alt == "C"))]
v[, read_strand := fifelse(ref == "A", "+", "-")]
nos <- fread(file.path(SP, "positions_noSNP.bed"), header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
v_nos <- v[nos, on = .(chr, pos), nomatch = NULL]
message("PGP1 SNVs ", n_vcf, "; PASS A>G / T>C ", nrow(v), "; at non-dbSNP-common candidate positions ", nrow(v_nos))

rd <- function(u) {
  fs <- list.files(file.path(SP, "pileups"), paste0("^", u, "__.*\\.tsv\\.gz$"), full.names = TRUE)
  d <- rbindlist(lapply(fs, function(f) fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", showProgress = FALSE)))
  d[v_nos[, .(chr, pos)], on = .(chr, pos), nomatch = NULL]
}
orient <- function(d) {
  x <- merge(v_nos[, .(chr, pos, ref, read_strand)], d[, !"ref"], by = c("chr", "pos"))
  x[, `:=`(alt_n = fifelse(read_strand == "+", G, c), ref_n = fifelse(read_strand == "+", ref_plus, ref_minus),
           third = fifelse(read_strand == "+", pmax(C, T), pmax(a, g)))]
  x[, cov := ref_n + alt_n]
  x[, .(chr, pos, cov, alt_n, third)]
}
ntc <- orient(rd("merged_NTC_ctrl"))
ko <- lapply(c("merged_ADAR1_KO", "merged_ADARB1_KO"), function(u) orient(rd(u))[, .(chr, pos, vaf = fifelse(cov >= 10, alt_n / cov, NA_real_))])
g <- merge(v_nos, ntc, by = c("chr", "pos"))
g[, multi := third >= MIN_THIRD_READS & third >= MIN_THIRD_FRAC * (cov + third)]
cnt <- data.table(step = c("PGP1 SNVs (all)", "PASS, A>G or T>C", "not dbSNP common, scanned candidate", "NTC reads on that orientation",
                           "cov >= 10", "alt >= 3", "VAF >= 5%", "not multiallelic"),
                  n = c(n_vcf, nrow(v), nrow(v_nos), nrow(g), g[cov >= MIN_COV, .N], g[cov >= MIN_COV & alt_n >= K, .N],
                        g[cov >= MIN_COV & alt_n >= K & alt_n / cov >= E, .N], g[cov >= MIN_COV & alt_n >= K & alt_n / cov >= E & !multi, .N]))
g <- g[cov >= MIN_COV & alt_n >= K & alt_n / cov >= E & !multi]
g <- merge(g, setnames(ko[[1]], "vaf", "ko1_vaf"), by = c("chr", "pos"), all.x = TRUE)
g <- merge(g, setnames(ko[[2]], "vaf", "ko2_vaf"), by = c("chr", "pos"), all.x = TRUE)
tr <- fread(file.path(PM, "gt_eval_v8", "truth", "gt_v7_truth_sites.tsv.gz"), select = "site")
g[, site := paste0(chr, ":", pos)]
stopifnot(!any(g$site %chin% tr$site))
out <- g[, .(set = "germline", site, chr, pos, genomic_ref = ref, genomic_alt = alt, read_strand, gt,
             ntc_cov = cov, ntc_alt = alt_n, ntc_vaf = round(alt_n / cov, 4), ko1_vaf = round(ko1_vaf, 4), ko2_vaf = round(ko2_vaf, 4))]
setorder(out, chr, pos)
fwrite(out, file.path(OUT, "germline_sites.tsv.gz"), sep = "\t")
cnt <- rbind(cnt, out[, .(step = paste0("kept, GT ", gt), n = .N), by = gt][, .(step, n)])
fwrite(cnt, file.path(OUT, "germline_counts.tsv"), sep = "\t")
print(cnt)
print(out[, .(n = .N, median_ntc_vaf = median(ntc_vaf), median_ko1_vaf = median(ko1_vaf, na.rm = TRUE),
              median_ko2_vaf = median(ko2_vaf, na.rm = TRUE)), by = gt])
