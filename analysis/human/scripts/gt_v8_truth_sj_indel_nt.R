#!/usr/bin/env Rscript
# gt_v8_truth_sj_indel_nt.R   (2026-10-24)
# KO ground truth sites (t = 50%, both gates; unique sites, independent of any caller) by distance (nt) to the nearest
# splice junction and to the nearest indel: positives vs informative negatives at 0, 1, ..., 10 and > 10 nt.
# Junctions and indels as in gt_v8_sj_indel_calls.R: STAR SJ.out.tab of the 8 NTC runs (>= 5 unique reads summed, or
# annotated; boundaries = first and last intron base); indels in the merged NTC reads with >= MIN_INDEL (3) reads
# (deletions: first and last deleted base; insertions: the reference bases on either side).
# Input : truth/gt_v7_truth_sites.tsv.gz, 01_starsolobams/*_SJ.out.tab, indel_sj/indels/<chr>.tsv.gz
# Output: human/pm/gt_eval_v8/indel_sj/truth_distance_nt.tsv  (feature x gate x class x d: n_sites)
#         human/pm/gt_eval_v8/indel_sj/truth_sites_distance.tsv.gz  (one row per gate x site)
#   Rscript human/scripts/gt_v8_truth_sj_indel_nt.R
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("SCRNA_DIR", ".")
PM <- file.path(B, "human/pm")
V8 <- file.path(PM, "gt_eval_v8")
OUT <- file.path(V8, "indel_sj")
MIN_INDEL <- as.integer(Sys.getenv("MIN_INDEL", "3"))
D_MAX <- 10L
D_LAB <- c(as.character(0:D_MAX), paste0("> ", D_MAX))

truth <- fread(file.path(V8, "truth/gt_v7_truth_sites.tsv.gz"), select = c("gate", "site", "chr", "pos", "class_t50"))
truth <- truth[class_t50 %chin% c("positive", "informative_negative")]

nearest <- function(pts, q) {
  setkey(pts, chr, b)
  qq <- q[, .(chr, b = pos, pos)]
  abs(pts[qq, on = .(chr, b), roll = "nearest", x.b] - qq$pos)
}
sjf <- list.files(file.path(PM, "01_starsolobams"), "_SJ\\.out\\.tab$", full.names = TRUE)
stopifnot(length(sjf) == 8)
sj <- rbindlist(lapply(sjf, fread, col.names = c("chr", "start", "end", "strand", "motif", "annot", "uniq", "multi",
                                                 "overhang")))
sj <- sj[, .(uniq = sum(uniq), annot = max(annot)), by = .(chr, start, end)][uniq >= 5 | annot == 1]
bd <- unique(rbind(sj[, .(chr, b = start)], sj[, .(chr, b = end)]))
inf <- list.files(file.path(OUT, "indels"), "\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(inf) == 23)
ind <- rbindlist(lapply(inf, fread, col.names = c("chr", "pos", "type", "len", "n_reads")))
ind <- ind[, .(n_reads = sum(n_reads)), by = .(chr, pos, type, len)][n_reads >= MIN_INDEL]
ib <- unique(rbind(ind[type == "D", .(chr, b = pos)], ind[type == "D", .(chr, b = pos + len - 1L)],
                   ind[type == "I", .(chr, b = pos - 1L)], ind[type == "I", .(chr, b = pos)]))

upos <- unique(truth[, .(chr, pos)])
upos[, `:=`(d_sj = nearest(bd, upos), d_indel = nearest(ib, upos))]
truth <- upos[truth, on = .(chr, pos)]
fwrite(truth, file.path(OUT, "truth_sites_distance.tsv.gz"), sep = "\t")

d_bin <- function(d) factor(fifelse(is.na(d) | d > D_MAX, paste0("> ", D_MAX), as.character(d)), D_LAB)
long <- rbind(truth[, .(feature = "junction", gate, class = class_t50, d = d_bin(d_sj))],
              truth[, .(feature = "indel", gate, class = class_t50, d = d_bin(d_indel))])
res <- long[, .(n_sites = .N), by = .(feature, gate, class, d)]
res <- res[CJ(feature = c("junction", "indel"), gate = c("expected", "depth"),
              class = c("positive", "informative_negative"), d = factor(D_LAB, D_LAB)),
           on = .(feature, gate, class, d)][is.na(n_sites), n_sites := 0L]
fwrite(res, file.path(OUT, "truth_distance_nt.tsv"), sep = "\t")
print(dcast(res, feature + gate + d ~ class, value.var = "n_sites")[order(feature, gate, d)])
