#!/usr/bin/env Rscript
# ============================================================================
# WASP reference-bias diagnostic — analysis + figures (Reviewer 2, Major 4)
# ============================================================================
# Compares strand-aware editing levels at REDIportal A-to-I sites between:
#   ALL reads (no WASP; standard STAR alignment)  = biased estimate
#   WASP-passed reads (vW==1 | no-vW)             = reference-bias-corrected
# If STAR reference bias underestimates editing, removing the allele-biased
# reads should RAISE editing levels (ref-supporting reads are preferentially
# removed). We quantify: (1) per-site editing-level shift, (2) allele
# composition of WASP-removed reads at editing sites (the direct bias
# signature), (3) coverage cost of WASP in single-cell data.
#
# Usage: Rscript wasp_refbias_analysis.R <OUTDIR> [MINCOV=20]
# Inputs (in OUTDIR): pileup_all.tsv, pileup_pass.tsv, readloss.tsv
# Outputs (in OUTDIR): wasp_refbias_summary.tsv, wasp_refbias_by_threshold.tsv,
#                      wasp_refbias_persite.tsv.gz, and PNG figures.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args   <- commandArgs(trailingOnly = TRUE)
OUT    <- if (length(args) >= 1) args[1] else "."
MINCOV <- if (length(args) >= 2) as.integer(args[2]) else 20L
theme_set(theme_bw(base_size = 12))

f_all  <- file.path(OUT, "pileup_all.tsv")
f_pass <- file.path(OUT, "pileup_pass.tsv")
stopifnot(file.exists(f_all), file.exists(f_pass))

read_pileup <- function(f) {
  d <- fread(f, sep = "\t", header = TRUE)
  # strand-aware columns from parse_pileup_stranded.pl:
  #   coverage = ref+alt at editing site; variant_cov = edited-allele count
  d[, editing_level := suppressWarnings(as.numeric(editing_level))]
  d[, ref_at_site := coverage - variant_cov]
  d[, .(chr, pos, ref, coverage, variant_cov, ref_at_site, editing_level)]
}

da <- read_pileup(f_all)
dp <- read_pileup(f_pass)
setnames(da, c("coverage","variant_cov","ref_at_site","editing_level"),
             c("cov_all","alt_all","ref_all","EL_all"))
setnames(dp, c("coverage","variant_cov","ref_at_site","editing_level"),
             c("cov_pass","alt_pass","ref_pass","EL_pass"))

# Base set: sites adequately covered WITHOUT WASP. Left-join WASP-passed; a site
# absent in the WASP pileup means all its reads were removed -> cov_pass = 0.
m <- merge(da, dp[, .(chr, pos, cov_pass, alt_pass, ref_pass, EL_pass)],
           by = c("chr","pos"), all.x = TRUE)
for (col in c("cov_pass","alt_pass","ref_pass"))
  m[is.na(get(col)), (col) := 0L]

fwrite(m, file.path(OUT, "wasp_refbias_persite.tsv.gz"), sep = "\t")

# ── Summaries across coverage thresholds ────────────────────────────────────
thresholds <- sort(unique(c(10L, 20L, 50L, MINCOV)))
summ_rows <- rbindlist(lapply(thresholds, function(th) {
  base <- m[cov_all >= th]                       # covered without WASP
  ev   <- base[cov_pass >= th & is.finite(EL_all) & is.finite(EL_pass)]  # paired-stable
  dEL  <- ev$EL_pass - ev$EL_all
  wp   <- tryCatch(wilcox.test(ev$EL_pass, ev$EL_all, paired = TRUE)$p.value,
                   error = function(e) NA_real_)
  data.table(
    min_cov            = th,
    n_sites_covered    = nrow(base),
    n_sites_evaluable  = nrow(ev),
    n_dropped_by_wasp  = sum(base$cov_pass < th),
    median_EL_all      = round(median(ev$EL_all), 4),
    median_EL_wasp     = round(median(ev$EL_pass), 4),
    median_dEL         = round(median(dEL), 4),
    mean_dEL           = round(mean(dEL), 4),
    pct_sites_EL_up    = round(100 * mean(dEL > 0), 2),
    pct_sites_EL_down  = round(100 * mean(dEL < 0), 2),
    wilcoxon_p         = wp
  )
}))
fwrite(summ_rows, file.path(OUT, "wasp_refbias_by_threshold.tsv"), sep = "\t")

# ── Allele composition of WASP-removed reads at editing sites (bias signature)
base <- m[cov_all >= MINCOV]
removed_ref <- pmax(base$ref_all - base$ref_pass, 0)
removed_alt <- pmax(base$alt_all - base$alt_pass, 0)
tot_rem_ref <- sum(removed_ref); tot_rem_alt <- sum(removed_alt)
tot_ref     <- sum(base$ref_all); tot_alt <- sum(base$alt_all)
# If reference bias exists, WASP removes ref-allele reads at a HIGHER rate than
# alt-allele reads at editing sites.
rate_rem_ref <- if (tot_ref > 0) tot_rem_ref / tot_ref else NA_real_
rate_rem_alt <- if (tot_alt > 0) tot_rem_alt / tot_alt else NA_real_

rl <- tryCatch(fread(file.path(OUT, "readloss.tsv")), error = function(e) NULL)

# ── Headline summary (long) ─────────────────────────────────────────────────
prim <- summ_rows[min_cov == MINCOV]
summary_long <- rbindlist(list(
  data.table(metric = "min_coverage",            value = as.character(MINCOV)),
  data.table(metric = "n_sites_evaluable",       value = as.character(prim$n_sites_evaluable)),
  data.table(metric = "median_EL_no_wasp",       value = as.character(prim$median_EL_all)),
  data.table(metric = "median_EL_wasp",          value = as.character(prim$median_EL_wasp)),
  data.table(metric = "median_delta_EL",         value = as.character(prim$median_dEL)),
  data.table(metric = "mean_delta_EL",           value = as.character(prim$mean_dEL)),
  data.table(metric = "pct_sites_EL_increased",  value = as.character(prim$pct_sites_EL_up)),
  data.table(metric = "pct_sites_EL_decreased",  value = as.character(prim$pct_sites_EL_down)),
  data.table(metric = "wilcoxon_paired_p",       value = format(prim$wilcoxon_p, scientific = TRUE)),
  data.table(metric = "editing_reads_removed_frac_at_sites", value = sprintf("%.4f", rate_rem_alt)),
  data.table(metric = "ref_reads_removed_frac_at_sites",     value = sprintf("%.4f", rate_rem_ref)),
  data.table(metric = "removed_ref_reads",       value = as.character(tot_rem_ref)),
  data.table(metric = "removed_alt_reads",       value = as.character(tot_rem_alt)),
  if (!is.null(rl)) data.table(metric = paste0("global_", rl$metric), value = as.character(rl$value))
), fill = TRUE)
fwrite(summary_long, file.path(OUT, "wasp_refbias_summary.tsv"), sep = "\t")

cat("\n==== WASP reference-bias summary (min_cov =", MINCOV, ") ====\n")
print(summary_long)
cat("\nBy threshold:\n"); print(summ_rows)

# ============================ FIGURES (ggplot2) =============================
ev <- m[cov_all >= MINCOV & cov_pass >= MINCOV & is.finite(EL_all) & is.finite(EL_pass)]
ev[, dEL := EL_pass - EL_all]
med_dEL <- median(ev$dEL)

# 1) EL_all vs EL_wasp hexbin with y=x
p1 <- ggplot(ev, aes(EL_all, EL_pass)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.5) +
  scale_fill_viridis_c(trans = "log10", name = "sites") +
  labs(title = "Editing level: no-WASP vs WASP-passed",
       subtitle = sprintf("REDIportal A-to-I sites, coverage >= %d in both (n = %s)",
                          MINCOV, format(nrow(ev), big.mark = ",")),
       x = "Editing level (all reads, no WASP)",
       y = "Editing level (WASP-passed reads)") +
  coord_equal(xlim = c(0,1), ylim = c(0,1))
ggsave(file.path(OUT, "wasp_refbias_scatter.png"), p1, width = 6, height = 5.5, dpi = 150)

# 2) Distribution of delta EL
p2 <- ggplot(ev, aes(dEL)) +
  geom_histogram(bins = 80, fill = "steelblue", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0, colour = "grey40", linetype = 2) +
  geom_vline(xintercept = med_dEL, colour = "red", linewidth = 0.6) +
  annotate("text", x = med_dEL, y = Inf, vjust = 1.5, hjust = -0.05,
           colour = "red", label = sprintf("median = %+.4f", med_dEL)) +
  labs(title = "Reference-bias correction per site",
       subtitle = "Positive = editing UNDER-estimated without WASP",
       x = expression(Delta*" editing level (WASP - noWASP)"), y = "sites")
ggsave(file.path(OUT, "wasp_refbias_deltaEL.png"), p2, width = 6.5, height = 4.5, dpi = 150)

# 3) Paired editing-level distributions
long <- rbind(
  data.table(EL = ev$EL_all,  grp = "no WASP"),
  data.table(EL = ev$EL_pass, grp = "WASP-passed"))
long[, grp := factor(grp, levels = c("no WASP","WASP-passed"))]
p3 <- ggplot(long, aes(grp, EL, fill = grp)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c("no WASP" = "grey60", "WASP-passed" = "steelblue"),
                    guide = "none") +
  labs(title = "Editing-level distribution", x = NULL, y = "Editing level")
ggsave(file.path(OUT, "wasp_refbias_violin.png"), p3, width = 5, height = 4.5, dpi = 150)

# 4) Allele composition of reads removed by WASP at editing sites
dbar <- data.table(
  allele = factor(c("reference (A/T)", "edited (G/C)"),
                  levels = c("reference (A/T)", "edited (G/C)")),
  removed = c(tot_rem_ref, tot_rem_alt),
  removal_rate = c(rate_rem_ref, rate_rem_alt))
p4 <- ggplot(dbar, aes(allele, removal_rate, fill = allele)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%\n(%s reads)", 100*removal_rate,
                                format(removed, big.mark = ","))),
            vjust = -0.2, size = 3.3) +
  scale_fill_manual(values = c("reference (A/T)" = "firebrick",
                               "edited (G/C)" = "seagreen"), guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "WASP-removed reads at editing sites, by allele",
       subtitle = "Higher reference removal = reference mapping bias",
       x = NULL, y = "fraction of allele's reads removed")
ggsave(file.path(OUT, "wasp_refbias_removed_alleles.png"), p4, width = 5.2, height = 4.5, dpi = 150)

# Combined panel if patchwork is available
ok <- requireNamespace("patchwork", quietly = TRUE)
if (ok) {
  library(patchwork)
  comb <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      title = "WASP reference-bias diagnostic — human PBMC (jzlxv3)",
      theme = theme(plot.title = element_text(face = "bold")))
  ggsave(file.path(OUT, "wasp_refbias_panel.png"), comb, width = 12, height = 10, dpi = 150)
}

cat("\nWrote figures + tables to", OUT, "\n")
