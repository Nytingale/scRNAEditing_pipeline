#!/usr/bin/env Rscript
# =============================================================================
# gt_truth_v7.R   (2026-10-08; A>G only from 2026-10-10)
#
# Teratoma ground truth v7: A-to-I truth positions.
#   truth position  chr:pos that clears the frozen v6 criteria (FROZEN_RULES.tsv: NTC editing,
#                   reproducibility across T1-T4, covered KO arm, KO depletion >= t) as RNA A>G
#                   in 10x read orientation (gt_truth_agtc_test_v6.R variant AG). Read-orientation
#                   T>C positions are not truth sites (user, 2026-10-10; the 2026-10-08 A>G + T>C
#                   version is in human/pm/_superseded/gt_eval_v7_agtc_2026-10-08/)
#   multiallelic    removed (third allele >= 3 reads and >= 5% of depth in merged NTC)
#   correct call    genomic A>G on the + strand at a genomic A, or genomic T>C on the - strand at
#                   a genomic T (= RNA-frame varType AG on either strand). Any other call there,
#                   including + strand T>C and - strand A>G, is wrong.
#   negatives       informative negatives of the same rows (confidently < 25% depleted in both
#                   KO arms); any call there is wrong.
# Gates: expected-gated (KO cov x NTC edit >= 3) and depth-gated (KO cov >= 10). t = 25/50/75%.
#
# Input : gt_eval_v6/agtc_test/truth_sites.tsv.gz (gt_truth_agtc_test_v6.R; job 1888116)
#         gt_eval_v6/agtc_test/truth_counts.tsv   (parity)
# Output: gt_eval_v7/truth/
#           gt_v7_truth_sites.tsv.gz   gate x eligible position: genomic_ref, a_to_i_strand,
#                                      read_type, read_strand, class_t25 / t50 / t75
#           gt_v7_truth_counts.tsv     class sizes per gate, with read_type split
#           gt_v7_<gate>_t<NN>_positives.bed
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
AG = Sys.getenv("GT_AGTC_DIR", file.path(PM, "gt_eval_v6/agtc_test"))
OUT = Sys.getenv("GT_V7_TRUTH_DIR", file.path(PM, "gt_eval_v7/truth")); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75)

x = fread(file.path(AG, "truth_sites.tsv.gz"))[variant == "AG" & eligible == TRUE & multiallelic == FALSE]
stopifnot(x[, uniqueN(site), by = gate][, all(V1 == x[, .N, by = gate]$N)], all(x$type == "AG"),
          all(x$strand %in% c("+", "-")))
# RNA A>G on + reads or T>C on - reads sits on a genomic A; the other two on a genomic T
x[, genomic_ref := fifelse((type == "AG") == (strand == "+"), "A", "T")]
x[, a_to_i_strand := fifelse(genomic_ref == "A", "+", "-")]
x[, c("chr", "pos") := tstrsplit(site, ":", fixed = TRUE)][, pos := as.integer(pos)]
# every truth site is an adenosine in the RNA: genomic A with the RNA on +, genomic T with it on -
reg = tempfile(fileext = ".txt"); u = unique(x[, .(chr, pos)])
writeLines(sprintf("%s:%d-%d", u$chr, u$pos, u$pos), reg)
fa = system2("samtools", c("faidx", "-r", reg, file.path(B, "human/ref_annots/hg38.analysisSet.fa")), stdout = TRUE)
u[, fasta_base := toupper(fa[!startsWith(fa, ">")])]
x = merge(x, u, by = c("chr", "pos"))
stopifnot(nrow(u) == sum(!startsWith(fa, ">")), all(x$fasta_base == x$genomic_ref))
message("FASTA check: ", nrow(u), " truth positions, all an adenosine on the RNA strand")
x[, fasta_base := NULL]
for (tn in names(DS)) {
  p = x[[paste0("pos_", tn)]]
  set(x, j = paste0("class_", tn), value = fifelse(p, "positive", fifelse(x$informative_negative, "informative_negative", "unresolved")))
}
stopifnot(!any(x$informative_negative & (x$pos_t25 | x$pos_t50 | x$pos_t75)))
tr = x[, c("gate", "site", "chr", "pos", "genomic_ref", "a_to_i_strand", "type", "strand", paste0("class_", names(DS))), with = FALSE]
setnames(tr, c("type", "strand"), c("read_type", "read_strand"))
setorder(tr, gate, chr, pos)
fwrite(tr, file.path(OUT, "gt_v7_truth_sites.tsv.gz"), sep = "\t")

cnt = rbindlist(lapply(names(DS), function(tn) {
  cl = paste0("class_", tn)
  tr[, .(t = tn, eligible = .N, positive = sum(get(cl) == "positive"),
         positive_read_AG = sum(get(cl) == "positive" & read_type == "AG"),
         positive_read_TC = sum(get(cl) == "positive" & read_type == "TC"),
         positive_genomic_A = sum(get(cl) == "positive" & genomic_ref == "A"),
         informative_negative = sum(get(cl) == "informative_negative"),
         unresolved = sum(get(cl) == "unresolved")), by = gate]
}))
ref = fread(file.path(AG, "truth_counts.tsv"))[variant == "AG" & multiallelic == "removed"]
chk = merge(cnt, ref[, .(gate, eligible_ref = eligible, informative_negative_ref = informative_negative,
                         t25 = pos_t25, t50 = pos_t50, t75 = pos_t75)], by = "gate")
chk[, positive_ref := fcase(t == "t25", t25, t == "t50", t50, t == "t75", t75)]
print(chk[, .(gate, t, eligible, eligible_ref, positive, positive_ref, informative_negative, informative_negative_ref)])
stopifnot(chk[, all(eligible == eligible_ref & positive == positive_ref & informative_negative == informative_negative_ref)])
fwrite(cnt, file.path(OUT, "gt_v7_truth_counts.tsv"), sep = "\t")
print(cnt)

for (g in unique(tr$gate)) for (tn in names(DS)) {
  b = tr[gate == g & get(paste0("class_", tn)) == "positive", .(chr, s = pos - 1L, e = pos, site, read_type, a_to_i_strand)]
  fwrite(b, file.path(OUT, sprintf("gt_v7_%s_%s_positives.bed", g, tn)), sep = "\t", col.names = FALSE)
}
cat("wrote", OUT, "\n")
