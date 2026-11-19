#!/usr/bin/env Rscript
# gt_v9_molecule_compare.R   (2026-11-17; v9 / 76-configuration copy of gt_v8_molecule_compare.R)
# Ranking with >= MIN_MOL variant molecules (molecule_filter_v8.R, scored as the headline: scoring_pwcov, gate_cov)
# against the headline ranking, per KO-based truth set; and the KO-based positives by variant molecules in the merged
# NTC arm (UMI-deduplicated BAM).
# Output: MW/compare_rank.tsv, MW/positives_molecules.tsv, MW/summary_callsets.tsv
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
V9 <- file.path(PM, "gt_eval_v9")
MW <- file.path(V9, "molecules")
COMP <- c(A = "T", C = "G", G = "C", T = "A")

main <- fread(file.path(V9, "scoring_pwcov76", "mismatch", "config_rank.tsv"))
mol <- fread(file.path(MW, "scoring_pwcov", "mismatch", "config_rank.tsv"))
x <- merge(main[, .(gate, t, config, rank_main = overall_rank, f1_main = mean_f1)],
           mol[, .(gate, t, config, rank_mol = overall_rank, f1_mol = mean_f1)], by = c("gate", "t", "config"))
stopifnot(x[, .N, by = .(gate, t)][, all(N == 76L)])
cmp <- x[, .(rho = cor(rank_main, rank_mol, method = "spearman"),
             top_main = config[which.min(rank_main)], top_main_rank_mol = rank_mol[which.min(rank_main)],
             top_mol = config[which.min(rank_mol)], median_d_f1 = median(f1_mol - f1_main),
             max_abs_d_f1 = max(abs(f1_mol - f1_main))), by = .(gate, t)]
fwrite(cmp, file.path(MW, "compare_rank.tsv"), sep = "\t")
print(cmp)

tr <- fread(file.path(V9, "truth", "gt_v9_truth_sites.tsv.gz"))
pos <- melt(tr, id.vars = c("gate", "site", "chr", "pos", "genomic_ref", "a_to_i_strand"),
            measure.vars = c("class_t25", "class_t50", "class_t75"), variable.name = "t", value.name = "class")[
  class == "positive"][, t := sub("class_", "", t)]
pos[, g_alt := ifelse(a_to_i_strand == "+", "G", "C")]
bams <- fread(file.path(MW, "bams.tsv"))
cn <- fread(bams[tier == "arm" & bam_type == "sampleUmiDup", counts])
p <- cn[pos, on = c("chr", "pos")]
stopifnot(!anyNA(p$ref))
p[, `:=`(alt_reads = fifelse(g_alt == "G", n_G, n_C), alt_mol = fifelse(g_alt == "G", m_G, m_C))]
pm <- p[, .(positives = .N, ge3_reads = sum(alt_reads >= 3), ge3_mol = sum(alt_mol >= 3),
            share_ge3_mol = mean(alt_mol >= 3), share_ge3_mol_of_ge3_reads = sum(alt_mol >= 3) / sum(alt_reads >= 3),
            reads_per_mol = sum(alt_reads) / sum(alt_mol)), by = .(gate, t)]
fwrite(pm, file.path(MW, "positives_molecules.tsv"), sep = "\t")
print(pm)

s <- rbindlist(lapply(c("arm", "sample", "celltype"), function(tr) fread(file.path(MW, "summary", paste0(tr, ".tsv")))[, tier := tr]),
               fill = TRUE)
ss <- s[n > 0, .(call_sets = .N, sites = sum(n), share_kept = sum(n_kept) / sum(n), share_ag_kept = sum(n_ag_kept) / sum(n_ag),
                 reads_per_mol_ag = sum(alt_reads_ag) / sum(alt_mol_ag), ref_mismatch = sum(ref_mismatch)),
        by = .(tier, caller, bam_type)]
fwrite(ss, file.path(MW, "summary_callsets.tsv"), sep = "\t")
print(ss[, .(share_ag_kept = sum(share_ag_kept * sites) / sum(sites)), by = tier])
