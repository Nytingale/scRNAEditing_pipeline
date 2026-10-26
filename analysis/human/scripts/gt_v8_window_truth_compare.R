#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_window_truth_compare.R   (2026-10-15; user: direct evidence on the clustered filter vs nomm)
# Both downstream_filter.R window passes applied to the ground-truth candidates, the way the truth
# construction applied the nearby-mismatch pass (gt_truth_v7_variants.R): all 12-type gate sites of a gate
# (read strand / type, PGP1 SNV positions removed) are the pool, and a filter is judged by what it keeps
# among the KO-depleted sites (d_max >= t, all 12 types; a real ADAR site is A-to-G, a depleted site of
# another type is an artifact).
#   none       all pool sites
#   mismatch   pass 1: removed if >= 2 mismatch types within w1 on its strand (truth's nomm_pass_pgp1, w1 = 200)
#   clustered  pass 2 only: kept if >= 2 pool sites (itself + >= 1 other, any type) within w2 on its strand
#   both       pass 1, then pass 2 among the pass-1 survivors (neighbours = survivors, so same type)
# w2 searched over 20..200 at each t by the highest A-to-G share of the depleted sites (the w1 rule, frac_weight = 1);
# the mismatch pass keeps the truth's w1 (asserted to also be the best at every t).
# Output: gt_eval_v8/truth/support/window_truth_compare.tsv   gate, t, state, w, n_depleted, n_ag, ag_pct,
#           ag_kept_pct (of none), n_neg, neg_kept_pct (A-to-G informative negatives, of none)
#         gt_eval_v8/truth/support/window_truth_grid.tsv      gate, t, pass, w, n_kept, n_ag, ag_pct, chosen (best w per t)
#         gt_eval_v8/truth/support/window_truth_compare_types.tsv  gate, t, state, type, n (KO-depleted sites kept)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
SUP <- file.path(PM, "gt_eval_v8", "truth", "support")
DS <- c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
W_ON <- seq(20, 200, 20)
W1 <- 200L

t12 <- fread(file.path(SUP, "truth12_site_flags.tsv.gz"))[pgp1 == FALSE]

# number of mismatch types among `nb` within +-w on the query's strand (gt_truth_v7_variants.R)
n_types_near <- function(q, nb, w) {
  q <- copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h <- nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL,
          .(row = i.row, type = x.type)]
  h <- rbind(h, q[, .(row, type)])
  h[, .(n_types = uniqueN(type)), by = row][order(row), n_types]
}

# number of pool sites (itself included) within +-w on the query's strand (filReliable sites_in_wind_n)
n_sites_near <- function(q, nb, w) {
  q <- copy(q)[, `:=`(row = .I, lo = pos - w, hi = pos + w)]
  h <- nb[q, on = .(strand, chr, pos >= lo, pos <= hi), allow.cartesian = TRUE, nomatch = NULL,
          .(row = i.row, site = paste(x.chr, x.pos))]
  n <- h[, .(n = uniqueN(site)), by = row]
  out <- integer(nrow(q)); out[n$row] <- n$n; out
}

res <- list(); grid <- list(); typ <- list()
for (g in c("expected", "depth")) {
  u <- t12[gate == g]
  stopifnot(nrow(u) > 0)
  keep <- list(none = rep(TRUE, nrow(u)), mismatch = u$nomm_pass_pgp1)
  wsel <- list(none = NA, mismatch = W1)
  for (ps in c("clustered", "both")) {
    base <- if (ps == "clustered") rep(TRUE, nrow(u)) else u$nomm_pass_pgp1
    pool <- u[base, .(chr, pos, strand)]
    near <- lapply(W_ON, function(w) { k <- base; k[base] <- n_sites_near(pool, pool, w) >= 2; k })
    gr <- rbindlist(lapply(names(DS), function(tn) rbindlist(lapply(seq_along(W_ON), function(i) {
      k <- near[[i]] & u$d_max >= DS[[tn]]
      data.table(gate = g, t = tn, pass = ps, w = W_ON[i], n_kept = sum(k), n_ag = sum(k & u$type == "AG"))
    }))))
    gr[, ag_pct := 100 * n_ag / pmax(n_kept, 1)]
    gr[, chosen := w == w[which.max(ag_pct)], by = t]
    grid[[paste(g, ps)]] <- gr
    keep[[ps]] <- lapply(setNames(names(DS), names(DS)), function(tn) near[[match(gr[t == tn & chosen == TRUE, w], W_ON)]])
    wsel[[ps]] <- lapply(setNames(names(DS), names(DS)), function(tn) gr[t == tn & chosen == TRUE, w])
  }
  # the mismatch pass's own w1 grid (truth construction), recomputed at every t; w1 chosen at t50
  m <- fread(file.path(SUP, "nomm_w1_grid.tsv"))[gate == g & pgp1 == TRUE]
  nb <- u[, .(chr, pos, strand, type)]
  mg <- rbindlist(lapply(names(DS), function(tn) {
    p <- u[d_max >= DS[[tn]], .(chr, pos, strand, type)]
    rbindlist(lapply(W_ON, function(w) {
      k <- n_types_near(p, nb, w) < 2
      data.table(gate = g, t = tn, pass = "mismatch", w = w, n_kept = sum(k), n_ag = sum(k & p$type == "AG"))
    }))
  }))
  stopifnot(identical(mg[t == "t50", n_kept], m$n_kept), identical(mg[t == "t50", n_ag], m$n_ag))
  mg[, ag_pct := 100 * n_ag / pmax(n_kept, 1)]
  mg[, chosen := w == w[which.max(ag_pct)], by = t]
  stopifnot(all(mg[chosen == TRUE, w] == m[chosen == TRUE, w1]))
  grid[[paste(g, "mismatch")]] <- mg
  for (st in names(keep)) for (tn in names(DS)) {
    dep <- u$d_max >= DS[[tn]]
    k <- if (is.list(keep[[st]])) keep[[st]][[tn]] else keep[[st]]
    neg <- u$type == "AG" & u$informative_negative
    res[[length(res) + 1]] <- data.table(
      gate = g, t = tn, state = st, w = if (is.list(wsel[[st]])) wsel[[st]][[tn]] else wsel[[st]], n_depleted = sum(k & dep), n_ag = sum(k & dep & u$type == "AG"),
      n_ag_none = sum(dep & u$type == "AG"), n_neg = sum(k & neg), n_neg_none = sum(neg))
    typ[[length(typ) + 1]] <- u[k & dep, .(n = .N), by = type][, `:=`(gate = g, t = tn, state = st)]
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
print(r[, .(gate, t, state, w, n_depleted, n_ag, ag_pct = round(ag_pct, 1), ag_kept_pct = round(ag_kept_pct, 1),
            n_neg, neg_kept_pct = round(neg_kept_pct, 1))])
