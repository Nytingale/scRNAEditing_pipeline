#!/usr/bin/env Rscript
# =============================================================================
# gt_truth_v9.R   (2026-11-07)
#
# Teratoma ground truth v9: KO-responsive A>G and T>C in read orientation, scored with
# A-to-G calls only (gt_score_v7.R GT_FP_RULE=ag).
#   positives   gt_truth_v7.R + gt_truth_v7_variants.R with GT_V7_VARIANT=AGTC: read-orientation A>G or T>C
#               positions clearing the frozen v6 criteria (counts within one read orientation; where both
#               orientations are candidates, the one with more merged-NTC reads), multiallelic removed,
#               cleaned as v7_nomm_pgp1 (12-type nearby-mismatch filter on the read strand, w1 searched on the
#               A>G share; PGP1 PASS + VQLOW SNV positions removed). The edited strand is the one where the
#               genomic base is A; a correct call is RNA-frame AG on that strand.
#   negatives   (same at every t) any AG call is a false positive
#     infneg      informative negatives of the same build before cleaning (confidently < 25% depleted in both
#                 KO arms; NTC editing and T1-T4 reproducibility as the positives)
#     pgp1_ag     PGP1 PASS SNV at a genomic A with a G allele or a genomic T with a C allele (any genotype):
#                 the G / C reads are the germline allele
#     pgp1_hom    PGP1 PASS SNV at a genomic A or T without the reference allele (1/1, 1/2, haploid 1) and
#                 not pgp1_ag: no adenosine on the edited strand. Heterozygous other-alt SNVs (0/1, ./1) are
#                 unscored (the reference allele can still be edited).
# Input : GT_V9_DIR/build/v7truth/gt_v7_truth_sites.tsv.gz (uncleaned), GT_V9_DIR/build/variants/v7_nomm_pgp1/,
#         pgp1_snps/PGP1_hu43860C_hg38_snvs.vcf.gz
# Output: GT_V9_DIR/truth/gt_v9_truth_sites.tsv.gz   gt_v7 columns + neg_source (gt_score_v7.R GT_V7_TRUTH)
#         GT_V9_DIR/truth/gt_v9_positive_sites.tsv.gz positive rows only (gt_v8_unit_cov.sh TRUTH)
#         GT_V9_DIR/truth/gt_v9_truth_counts.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
V9 = Sys.getenv("GT_V9_DIR", file.path(PM, "gt_eval_v9"))
OUT = file.path(V9, "truth"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DS = c("t25", "t50", "t75"); CL = paste0("class_", DS)
COLS = c("gate", "site", "chr", "pos", "genomic_ref", "a_to_i_strand", "read_type", "read_strand", CL)

cln = fread(file.path(V9, "build/variants/v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz"))
unc = fread(file.path(V9, "build/v7truth/gt_v7_truth_sites.tsv.gz"))
stopifnot(identical(names(cln), COLS), identical(names(unc), COLS), all(cln$read_type %in% c("AG", "TC")))
GATES = unique(cln$gate)

# ---- negatives -----------------------------------------------------------------------------------
infneg = unc[class_t25 == "informative_negative", .(gate, site, chr, pos, genomic_ref, neg_source = "infneg")]
stopifnot(unc[class_t25 == "informative_negative", all(class_t50 == class_t25 & class_t75 == class_t25)])

vcf = fread(cmd = paste("gzip -cd", shQuote(file.path(PM, "pgp1_snps/PGP1_hu43860C_hg38_snvs.vcf.gz")), "| grep -v '^#'"),
            header = FALSE, select = c(1, 2, 4, 5, 7, 10), col.names = c("chr", "pos", "ref", "alt", "filter", "gt"))
vcf = vcf[filter == "PASS" & ref %in% c("A", "T")]
vcf[, gt := sub(":.*", "", gt)]
al = strsplit(vcf$alt, ",", fixed = TRUE)
gtal = strsplit(vcf$gt, "[/|]")
called = mapply(function(a, g) { i = suppressWarnings(as.integer(g)); a[i[!is.na(i) & i > 0]] }, al, gtal, SIMPLIFY = FALSE)
no_ref = vapply(gtal, function(g) all(g != "0" & g != "."), logical(1))
edit_base = c(A = "G", T = "C")
vcf[, has_edit_allele := mapply(function(b, r) edit_base[[r]] %in% b, called, ref)]
vcf[, src := fifelse(has_edit_allele, "pgp1_ag", fifelse(no_ref, "pgp1_hom", NA_character_))]
print(vcf[, .N, by = .(gt, src)][order(src, -N)])
pg = vcf[!is.na(src), .(site = paste0(chr, ":", pos), chr, pos, genomic_ref = ref, neg_source = src)]
neg = rbind(infneg, rbindlist(lapply(GATES, function(g) copy(pg)[, gate := g])), use.names = TRUE)
# a position in more than one negative group keeps all labels
neg = neg[, .(chr = chr[1], pos = pos[1], genomic_ref = genomic_ref[1],
              neg_source = paste(sort(unique(neg_source)), collapse = "+")), by = .(gate, site)]
stopifnot(neg[, uniqueN(genomic_ref), by = .(gate, site)][, all(V1 == 1)])

# negatives must not be positives at any t (PGP1 positions are removed from the positives by the cleaning)
pos_sites = cln[class_t25 == "positive" | class_t50 == "positive" | class_t75 == "positive", .(gate, site)]
clash = merge(neg, pos_sites, by = c("gate", "site"))
stopifnot(nrow(clash) == 0L)

# ---- assemble ------------------------------------------------------------------------------------
# cleaned rows keep their classes except informative negatives, which come from `neg`; rows of `neg` replace
# any cleaned row at the same gate x site
keep = cln[!(class_t25 == "informative_negative")][, neg_source := NA_character_]
keep = keep[!neg, on = c("gate", "site")]
negr = neg[, .(gate, site, chr, pos, genomic_ref, a_to_i_strand = fifelse(genomic_ref == "A", "+", "-"),
               read_type = NA_character_, read_strand = NA_character_,
               class_t25 = "informative_negative", class_t50 = "informative_negative", class_t75 = "informative_negative",
               neg_source)]
tr = rbind(keep, negr, use.names = TRUE)
setorder(tr, gate, chr, pos)
stopifnot(tr[, .N, by = .(gate, site)][, all(N == 1L)])
fwrite(tr, file.path(OUT, "gt_v9_truth_sites.tsv.gz"), sep = "\t")
fwrite(tr[class_t25 == "positive" | class_t50 == "positive" | class_t75 == "positive"],
       file.path(OUT, "gt_v9_positive_sites.tsv.gz"), sep = "\t")

cnt = rbindlist(lapply(DS, function(tn) {
  cl = paste0("class_", tn)
  tr[, .(t = tn, positive = sum(get(cl) == "positive"),
         positive_read_AG = sum(get(cl) == "positive" & read_type %in% "AG"),
         positive_read_TC = sum(get(cl) == "positive" & read_type %in% "TC"),
         negative = sum(get(cl) == "informative_negative"),
         neg_infneg = sum(get(cl) == "informative_negative" & grepl("infneg", neg_source)),
         neg_pgp1_ag = sum(get(cl) == "informative_negative" & grepl("pgp1_ag", neg_source)),
         neg_pgp1_hom = sum(get(cl) == "informative_negative" & grepl("pgp1_hom", neg_source)),
         unresolved = sum(get(cl) == "unresolved")), by = gate]
}))
fwrite(cnt, file.path(OUT, "gt_v9_truth_counts.tsv"), sep = "\t")
print(cnt)
cat("wrote", OUT, "\n")
