#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_window_cluster_retention.R   (2026-10-18; R2 major 3: does the nearby-mismatch pass discard hyper-edited clusters?)
# Pool = truth12 gate sites, PGP1 SNV positions removed (as gt_v8_window_truth_compare.R). For the KO-depleted
# A-to-G sites (d_max >= t), count the other A-to-G pool sites within +-200 nt on the same strand (w1 = 200), bin
# 0 / 1 / 2 / 3 / > 3, and report the share kept by pass 1 (nomm_pass_pgp1: < 2 mismatch types within w1).
# Output: gt_eval_v8/truth/support/window_cluster_retention.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
SUP <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8/truth/support")
DS <- c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
W1 <- 200L

t12 <- fread(file.path(SUP, "truth12_site_flags.tsv.gz"))[pgp1 == FALSE]

n_ag_near <- function(q, nb, w) {
  q <- copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h <- nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL,
          .(row = i.row, p = x.pos, qp = i.pos)][p != qp]
  out <- integer(nrow(q)); n <- h[, .N, by = row]; out[n$row] <- n$N; out
}

res <- rbindlist(lapply(c("expected", "depth"), function(g) {
  u <- t12[gate == g]
  ag <- u[type == "AG", .(chr, pos, strand)]
  rbindlist(lapply(names(DS), function(tn) {
    d <- u[type == "AG" & d_max >= DS[[tn]]]
    d[, n_near := n_ag_near(d[, .(chr, pos, strand)], ag, W1)]
    d[, cl := fifelse(n_near > 3, ">3", as.character(n_near))]
    rbind(d[, .(n = .N, kept = sum(nomm_pass_pgp1)), by = cl],
          d[, .(cl = "all", n = .N, kept = sum(nomm_pass_pgp1))])[, `:=`(gate = g, t = tn)]
  }))
}))
res[, kept_pct := round(100 * kept / n, 1)]
setcolorder(res, c("gate", "t", "cl"))
fwrite(res, file.path(SUP, "window_cluster_retention.tsv"), sep = "\t")
print(res[order(gate, t, match(cl, c("0", "1", "2", "3", ">3", "all")))])
