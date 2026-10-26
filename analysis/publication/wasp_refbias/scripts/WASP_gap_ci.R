#!/usr/bin/env Rscript
# ============================================================================
# WASP_gap_ci.R — how uncertain IS the reference-vs-edited removal gap?
# ============================================================================
# >>> RUN VIA WASP_gap_ci.sh (SGE). Reads ~1 GB of pileups and bootstraps over
# >>> ~1.8 M sites. NOT a login-node script.
#
# WHY THIS EXISTS. WASP_prep_source_data.R reports the gap's CI from
# prop.test(), which treats every READ as an independent observation. Reads are
# not independent: reads at one site share the same local sequence, the same
# REDIportal variant density and the same mapping context, so whether WASP
# removes them is strongly correlated within a site. With ~300 M reference
# reads the nominal CI is +/-0.03 pp, which is not credible, and any reviewer
# who asks "is that difference significant?" is entitled to reject the answer.
#
# THE SITE IS THE UNIT. Three intervals are computed per threshold so they can
# be compared directly:
#   naive       prop.test over reads          (what stage 2 currently reports)
#   linearised  cluster-robust, sites         (one pass, closed form)
#   bootstrap   resample sites with replacement (independent check on the above)
#
# Linearisation for the difference of two ratio estimators sharing clusters.
# Per site i: a_i/b_i = ref removed/total, c_i/d_i = edited removed/total.
#   p_ref = sum(a)/sum(b),  p_alt = sum(c)/sum(d),  theta = p_ref - p_alt
#   u_i   = (a_i - p_ref*b_i)/sum(b) - (c_i - p_alt*d_i)/sum(d)
#   SE    = sqrt( n/(n-1) * sum(u_i^2) )        [sum(u_i) = 0 by construction]
#
# Writes: report/<dataset>/WASP_gap_ci.tsv
#
# Usage: Rscript WASP_gap_ci.R [dataset ...] [--reps N]
# ============================================================================
suppressPackageStartupMessages(library(data.table))

this_dir <- local({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
})
DELIV_DIR <- normalizePath(file.path(this_dir, ".."), mustWork = FALSE)
REPO_DIR  <- normalizePath(file.path(DELIV_DIR, "..", ".."), mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)
# 1000 replicates. NOTE this is the number of RESAMPLES, not a subsample of
# sites — every replicate draws the full site set (1.77-1.88 M) with
# replacement. Cost is linear in reps: 500 took 558 s for 2 datasets x 3
# thresholds including ~3 min of pileup reading, at 3.33 GB peak. Reps only
# stabilise the percentile endpoints; the SE is already confirmed by the
# closed-form linearisation, which agrees to within 1-6%.
reps <- 1000L
if (any(args == "--reps")) {
  i <- which(args == "--reps")
  reps <- as.integer(args[i + 1L]); args <- args[-c(i, i + 1L)]
}
datasets <- if (length(args)) args else c("jzlxv3", "jzlss2")

reg <- fread(file.path(this_dir, "WASP_datasets.tsv"))
set.seed(8787L)   # manuscript_seed

read_pileup <- function(f) {
  d <- fread(f, sep = "\t", header = TRUE,
             select = c("chr", "pos", "coverage", "variant_cov"))
  d[, ref_at_site := coverage - variant_cov]
  d[]
}

for (ds in datasets) {
  cfg <- reg[dataset == ds]
  stopifnot(nrow(cfg) == 1L)
  A <- file.path(REPO_DIR, "human", ds, "wasp_refbias")
  cat("===", ds, "(", cfg$label, ") ===\n"); flush.console()

  da <- read_pileup(file.path(A, "pileup_all.tsv"))
  dp <- read_pileup(file.path(A, "pileup_pass.tsv"))
  setnames(da, c("coverage","variant_cov","ref_at_site"), c("cov_all","alt_all","ref_all"))
  setnames(dp, c("coverage","variant_cov","ref_at_site"), c("cov_pass","alt_pass","ref_pass"))
  m <- merge(da, dp[, .(chr, pos, cov_pass, alt_pass, ref_pass)],
             by = c("chr", "pos"), all.x = TRUE)
  for (col in c("cov_pass","alt_pass","ref_pass")) m[is.na(get(col)), (col) := 0L]
  rm(da, dp); invisible(gc())

  out <- rbindlist(lapply(c(10L, 20L, 50L), function(th) {
    b  <- m[cov_all >= th]
    a  <- as.numeric(pmax(b$ref_all - b$ref_pass, 0)); bb <- as.numeric(b$ref_all)
    cc <- as.numeric(pmax(b$alt_all - b$alt_pass, 0)); dd <- as.numeric(b$alt_all)
    n  <- length(a)
    Sa <- sum(a); Sb <- sum(bb); Sc <- sum(cc); Sd <- sum(dd)
    p_ref <- Sa / Sb; p_alt <- Sc / Sd
    theta <- p_ref - p_alt

    se_naive <- sqrt(p_ref * (1 - p_ref) / Sb + p_alt * (1 - p_alt) / Sd)
    u <- (a - p_ref * bb) / Sb - (cc - p_alt * dd) / Sd
    se_lin <- sqrt(n / (n - 1) * sum(u^2))

    th_b <- numeric(reps)
    for (r in seq_len(reps)) {
      i <- sample.int(n, n, replace = TRUE)
      th_b[r] <- sum(a[i]) / sum(bb[i]) - sum(cc[i]) / sum(dd[i])
    }
    qb <- unname(quantile(th_b, c(0.025, 0.975)))
    z  <- theta / se_lin

    cat(sprintf(
      "  cov>=%-2d gap=%6.3f pp | naive SE %5.3f | linearised SE %5.3f CI[%6.3f,%6.3f] p=%-8s | boot SE %5.3f CI[%6.3f,%6.3f] | inflation %.1fx\n",
      th, 100*theta, 100*se_naive, 100*se_lin,
      100*(theta - 1.96*se_lin), 100*(theta + 1.96*se_lin),
      format.pval(2 * pnorm(-abs(z)), digits = 3),
      100*sd(th_b), 100*qb[1], 100*qb[2], se_lin / se_naive)); flush.console()

    data.table(
      dataset = ds, min_cov = th, n_sites = n,
      reads_ref = Sb, reads_alt = Sd,
      gap_pp = round(100 * theta, 4),
      se_naive_pp = round(100 * se_naive, 4),
      naive_lo_pp = round(100 * (theta - 1.96 * se_naive), 4),
      naive_hi_pp = round(100 * (theta + 1.96 * se_naive), 4),
      se_cluster_pp = round(100 * se_lin, 4),
      cluster_lo_pp = round(100 * (theta - 1.96 * se_lin), 4),
      cluster_hi_pp = round(100 * (theta + 1.96 * se_lin), 4),
      cluster_p = 2 * pnorm(-abs(z)),
      boot_reps = reps,
      se_boot_pp = round(100 * sd(th_b), 4),
      boot_lo_pp = round(100 * qb[1], 4),
      boot_hi_pp = round(100 * qb[2], 4),
      boot_frac_le0 = round(mean(th_b <= 0), 4),
      se_inflation = round(se_lin / se_naive, 2)
    )
  }))

  rdir <- file.path(DELIV_DIR, "report", ds)
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(out, file.path(rdir, "WASP_gap_ci.tsv"), sep = "\t")
  cat("  wrote", file.path("report", ds, "WASP_gap_ci.tsv"), "\n\n"); flush.console()
  rm(m); invisible(gc())
}
