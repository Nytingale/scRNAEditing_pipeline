#!/usr/bin/env Rscript
# =============================================================================
# gt_v9_window_cluster_compare.R   (2026-11-16; Response Figure 12 on the v9 truth, teratoma only: the aorta is held out)
# v9 port of gt_v8_window_truth_compare.R + gt_v8_window_cluster_retention.R. Pool = all 12-type gate sites of a gate
# (read strand / type, PGP1 SNV positions removed; gt_eval_v9/build/variants/truth12_site_flags.tsv.gz). A filter is
# judged on the KO-depleted sites (d_max >= t, all 12 types). v9: A-to-I = A>G or T>C in read orientation (both are
# positives); columns keep the v8 names, so n_ag / ag_pct / ag_kept_pct count A-to-I.
#   none       all pool sites
#   mismatch   pass 1: removed if >= 2 mismatch types within w1 on its strand (truth's nomm_pass_pgp1; w1 per gate
#              from nomm_w1_grid.tsv, chosen at t50 on the A>G share by the truth build)
#   clustered  pass 2 only: kept if >= 2 pool sites (itself + >= 1 other, any type) within w2 on its strand
#   both       pass 1, then pass 2 among the pass-1 survivors
# w2 searched over 20..200 at each t by the highest A-to-I share of the depleted sites.
# Clusters: KO-depleted A-to-I sites binned by other A-to-I pool sites within w1 on the same strand
# (0 / 1 / 2 / 3 / > 3), share kept by pass 1.
# Output: gt_eval_v9/truth/support/{window_truth_compare.tsv, window_truth_grid.tsv, window_truth_compare_types.tsv,
#         window_cluster_retention.tsv}  (formats as v8; read by GTDEF_v8_nomm_cluster_retention_source_data.R)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
V <- file.path(PM, "gt_eval_v9", "build", "variants")
SUP <- file.path(PM, "gt_eval_v9", "truth", "support")
dir.create(SUP, recursive = TRUE, showWarnings = FALSE)
DS <- c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
W_ON <- seq(20, 200, 20)
AI <- c("AG", "TC")
W1GRID <- fread(file.path(V, "nomm_w1_grid.tsv"))
W1G <- W1GRID[pgp1 == TRUE & chosen == TRUE, setNames(w1, gate)]

t12 <- fread(file.path(V, "truth12_site_flags.tsv.gz"))[pgp1 == FALSE]
t12[, ai := type %in% AI]

n_types_near <- function(q, nb, w) {
  q <- copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h <- nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL,
          .(row = i.row, type = x.type)]
  h <- rbind(h, q[, .(row, type)])
  h[, .(n_types = uniqueN(type)), by = row][order(row), n_types]
}
n_sites_near <- function(q, nb, w) {
  q <- copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h <- nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL,
          .(row = i.row, site = paste(x.chr, x.pos))]
  n <- h[, .(n = uniqueN(site)), by = row]
  out <- integer(nrow(q)); out[n$row] <- n$n; out
}
n_other_near <- function(q, nb, w) {
  q <- copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h <- nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL,
          .(row = i.row, p = x.pos, qp = i.pos)][p != qp]
  out <- integer(nrow(q)); n <- h[, .N, by = row]; out[n$row] <- n$N; out
}

res <- list(); grid <- list(); typ <- list(); clus <- list()
for (g in c("expected", "depth")) {
  u <- t12[gate == g]
  stopifnot(nrow(u) > 0)
  nb <- u[, .(chr, pos, strand, type)]
  stopifnot(identical(n_types_near(nb, nb, W1G[[g]]) < 2, u$nomm_pass_pgp1))
  keep <- list(none = rep(TRUE, nrow(u)), mismatch = u$nomm_pass_pgp1)
  wsel <- list(none = NA, mismatch = W1G[[g]])
  for (ps in c("clustered", "both")) {
    base <- if (ps == "clustered") rep(TRUE, nrow(u)) else u$nomm_pass_pgp1
    pool <- u[base, .(chr, pos, strand)]
    near <- lapply(W_ON, function(w) { k <- base; k[base] <- n_sites_near(pool, pool, w) >= 2; k })
    gr <- rbindlist(lapply(names(DS), function(tn) rbindlist(lapply(seq_along(W_ON), function(i) {
      k <- near[[i]] & u$d_max >= DS[[tn]]
      data.table(gate = g, t = tn, pass = ps, w = W_ON[i], n_kept = sum(k), n_ag = sum(k & u$ai))
    }))))
    gr[, ag_pct := 100 * n_ag / pmax(n_kept, 1)]
    gr[, chosen := w == w[which.max(ag_pct)], by = t]
    grid[[paste(g, ps)]] <- gr
    keep[[ps]] <- lapply(setNames(names(DS), names(DS)), function(tn) near[[match(gr[t == tn & chosen == TRUE, w], W_ON)]])
    wsel[[ps]] <- lapply(setNames(names(DS), names(DS)), function(tn) gr[t == tn & chosen == TRUE, w])
  }
  # the mismatch pass's own w grid at every t (A-to-I share); parity with the truth build's grid at t50 (A>G counts)
  m <- W1GRID[gate == g & pgp1 == TRUE]
  mg <- rbindlist(lapply(names(DS), function(tn) {
    p <- u[d_max >= DS[[tn]], .(chr, pos, strand, type)]
    rbindlist(lapply(W_ON, function(w) {
      k <- n_types_near(p, nb, w) < 2
      data.table(gate = g, t = tn, pass = "mismatch", w = w, n_kept = sum(k), n_ag = sum(k & p$type %in% AI),
                 n_ag_only = sum(k & p$type == "AG"))
    }))
  }))
  stopifnot(identical(mg[t == "t50", n_kept], m$n_kept), identical(mg[t == "t50", n_ag_only], m$n_ag))
  mg[, n_ag_only := NULL]
  mg[, ag_pct := 100 * n_ag / pmax(n_kept, 1)]
  mg[, chosen := w == w[which.max(ag_pct)], by = t]
  grid[[paste(g, "mismatch")]] <- mg
  neg <- u$ai & u$informative_negative
  for (st in names(keep)) for (tn in names(DS)) {
    dep <- u$d_max >= DS[[tn]]
    k <- if (is.list(keep[[st]])) keep[[st]][[tn]] else keep[[st]]
    res[[length(res) + 1]] <- data.table(
      gate = g, t = tn, state = st, w = if (is.list(wsel[[st]])) wsel[[st]][[tn]] else wsel[[st]],
      n_depleted = sum(k & dep), n_ag = sum(k & dep & u$ai), n_ag_none = sum(dep & u$ai),
      n_neg = sum(k & neg), n_neg_none = sum(neg))
    typ[[length(typ) + 1]] <- u[k & dep, .(n = .N), by = type][, `:=`(gate = g, t = tn, state = st)]
  }
  aip <- u[ai == TRUE, .(chr, pos, strand)]
  for (tn in names(DS)) {
    d <- u[ai == TRUE & d_max >= DS[[tn]]]
    d[, n_near := n_other_near(d[, .(chr, pos, strand)], aip, W1G[[g]])]
    d[, cl := fifelse(n_near > 3, ">3", as.character(n_near))]
    clus[[length(clus) + 1]] <- rbind(d[, .(n = .N, kept = sum(nomm_pass_pgp1)), by = cl],
                                      d[, .(cl = "all", n = .N, kept = sum(nomm_pass_pgp1))])[, `:=`(gate = g, t = tn)]
  }
}
r <- rbindlist(res)
r[, `:=`(ag_pct = 100 * n_ag / n_depleted, ag_kept_pct = 100 * n_ag / n_ag_none, neg_kept_pct = 100 * n_neg / n_neg_none)]
fwrite(r, file.path(SUP, "window_truth_compare.tsv"), sep = "\t")
fwrite(rbindlist(grid, use.names = TRUE), file.path(SUP, "window_truth_grid.tsv"), sep = "\t")
ty <- rbindlist(typ, use.names = TRUE)[, .(gate, t, state, type, n)]
stopifnot(identical(ty[, .(n_depleted = sum(n)), by = .(gate, t, state)][r[, .(gate, t, state)], on = .(gate, t, state), n_depleted],
                    r$n_depleted))
fwrite(ty, file.path(SUP, "window_truth_compare_types.tsv"), sep = "\t")
cr <- rbindlist(clus)[, kept_pct := round(100 * kept / n, 1)]
setcolorder(cr, c("gate", "t", "cl"))
fwrite(cr, file.path(SUP, "window_cluster_retention.tsv"), sep = "\t")
r[, oth_kept_pct := 100 * (n_depleted - n_ag) / (n_depleted[state == "none"] - n_ag[state == "none"]), by = .(gate, t)]
print(r[, .(gate, t, state, w, n_depleted, n_ag, ag_pct = round(ag_pct, 1), ag_kept_pct = round(ag_kept_pct, 1),
            oth_kept_pct = round(oth_kept_pct, 1), n_neg, neg_kept_pct = round(neg_kept_pct, 1))])
print(cr[order(gate, t, match(cl, c("0", "1", "2", "3", ">3", "all")))])
