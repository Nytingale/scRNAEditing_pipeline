#!/usr/bin/env Rscript
# =============================================================================
# gt_eval_tier_errorbars.R
#
# Uncertainty for the teratoma GT method rankings, per tier:
#
#   sample   biological SE across the four teratomas (T1–T4). Discovery is
#            restricted to arms of one teratoma at a time; each is scored
#            against the FIXED v3 truth set. Mean ± SE (n = 4).
#
#   dataset  no biological replicate axis (one merged BAM per arm). Site-level
#            bootstrap SE: resample the combo's discovered sites with
#            replacement, recompute precision / recall / F1 (B = 500).
#
# Writes (under --outdir, default true_def_v4_full):
#   gt_tier_scoring_with_se.tsv   pooled estimate + SE columns
#   gt_tier_replicate_scores.tsv  per-teratoma scores (sample tier only)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
parse_cli <- function() {
  a <- commandArgs(trailingOnly = TRUE); out <- list(); i <- 1L
  while (i <= length(a)) {
    if (grepl("^--", a[i])) { out[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L
    } else i <- i + 1L
  }
  out
}
args <- parse_cli()

B        <- as.integer(args[["bootstrap"]] %||% 500L)
SEED     <- as.integer(args[["seed"]]      %||% 20260919L)
outdir   <- args[["outdir"]] %||%
  file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval/true_def_v4_full")
truth_f  <- args[["truth"]] %||% file.path(outdir, "gt_truth_sites.tsv.gz")
disc_spec <- args[["discovery"]] %||% paste0(
  paste0("sample=", Sys.getenv("SCRNA_DIR", "."), "/human/pm/gt_eval/persample_union_v4/gt_eval_discovery.tsv.gz,"),
  paste0("dataset=", Sys.getenv("SCRNA_DIR", "."), "/human/pm/gt_eval/merged_union/gt_eval_discovery_balanced.tsv.gz,"),
  paste0("celltype=", Sys.getenv("SCRNA_DIR", "."), "/human/pm/gt_eval/celltype_union/gt_eval_discovery_poolB.tsv.gz")
)
fp_classes <- c("FP_novel", "FP_known_not_depleted")
# truth_scope: union = global TRUE count (v3/v4); per_tier = TRUE ∩ this
# tier's discovery sites (v5). Replicate / bootstrap SEs use the same
# yardstick as the pooled score for that tier.
truth_scope <- args[["truth_scope"]] %||% "union"
if (!truth_scope %in% c("union", "per_tier"))
  stop("--truth_scope must be 'union' or 'per_tier'")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

score_sites <- function(sites, truth_key, n_true_total) {
  if (!length(sites)) {
    return(list(n_discovered = 0L, n_true = 0L, n_fp = 0L, n_evaluable = 0L,
                precision = NA_real_, recall_vs_union = NA_real_, f1 = NA_real_))
  }
  dk <- unique(data.table(site = sites))
  dk <- merge(dk, truth_key, by = "site", all.x = TRUE)
  dk[is.na(class), class := "not_in_truth"]
  n_tr <- sum(dk$class == "TRUE")
  n_fp <- sum(dk$class %in% fp_classes)
  n_eval <- n_tr + n_fp
  prec <- if (n_eval > 0) n_tr / n_eval else NA_real_
  rec  <- if (n_true_total > 0) n_tr / n_true_total else NA_real_
  f1   <- if (is.finite(prec) && is.finite(rec) && (prec + rec) > 0)
            2 * prec * rec / (prec + rec) else NA_real_
  list(n_discovered = nrow(dk), n_true = n_tr, n_fp = n_fp, n_evaluable = n_eval,
       precision = prec, recall_vs_union = rec, f1 = f1)
}

msg("loading truth: ", truth_f)
truth <- fread(cmd = paste("gzip -cd", shQuote(truth_f)), sep = "\t",
               showProgress = FALSE)
truth_key <- truth[, .(site, class)]
n_true_global <- sum(truth$class == "TRUE")
msg("TRUE sites in global yardstick: ", n_true_global,
    "  truth_scope=", truth_scope)

pooled <- fread(file.path(outdir, "gt_tier_scoring.tsv"))
pooled[, f1 := fifelse(precision + recall_vs_union > 0,
                       2 * precision * recall_vs_union / (precision + recall_vs_union),
                       NA_real_)]

rep_rows <- list()
se_rows  <- list()

specs <- strsplit(strsplit(disc_spec, ",")[[1]], "=")
set.seed(SEED)

for (s in specs) {
  if (length(s) != 2L) next
  tier_name <- s[1]; path <- s[2]
  if (!file.exists(path)) { msg("SKIP ", tier_name, ": missing ", path); next }
  msg("tier ", tier_name, ": reading ", path)
  d <- fread(cmd = paste("gzip -cd", shQuote(path)), sep = "\t",
             showProgress = FALSE)
  gdims <- intersect(c("caller", "bam_type", "strand_mode", "snp_subtract"),
                     names(d))
  msg("  ", nrow(d), " rows, ", uniqueN(d[, gdims, with = FALSE]), " combos")

  if (truth_scope == "per_tier") {
    n_true_total <- sum(truth$site %in% unique(d$site) & truth$class == "TRUE")
    msg("  per-tier n_true_total=", n_true_total,
        " (of ", n_true_global, " global)")
  } else {
    n_true_total <- n_true_global
  }

  if (tier_name == "sample" && "arm" %in% names(d)) {
    d[, teratoma := sub("_.*$", "", arm)]
    tier_reps <- list()
    for (tx in sort(unique(d$teratoma))) {
      msg("  teratoma ", tx)
      dx <- unique(d[teratoma == tx, c(gdims, "site"), with = FALSE])
      sc <- dx[, {
        r <- score_sites(site, truth_key, n_true_total)
        as.data.table(r)
      }, by = gdims]
      sc[, `:=`(tier = tier_name, teratoma = tx)]
      tier_reps[[length(tier_reps) + 1L]] <- sc
    }
    reps <- rbindlist(tier_reps, fill = TRUE)
    rep_rows[[length(rep_rows) + 1L]] <- reps
    sum_se <- reps[, {
      .(n_reps = .N,
        precision_mean = mean(precision, na.rm = TRUE),
        precision_se   = sd(precision,   na.rm = TRUE) / sqrt(.N),
        recall_mean    = mean(recall_vs_union, na.rm = TRUE),
        recall_se      = sd(recall_vs_union,   na.rm = TRUE) / sqrt(.N),
        f1_mean        = mean(f1, na.rm = TRUE),
        f1_se          = sd(f1,   na.rm = TRUE) / sqrt(.N))
    }, by = gdims]
    sum_se[, `:=`(tier = tier_name, se_method = "teratoma_replicate")]
    se_rows[[length(se_rows) + 1L]] <- sum_se

  } else if (tier_name == "dataset") {
    msg("  bootstrap B=", B)
    dk <- unique(d[, c(gdims, "site"), with = FALSE])
    dk <- merge(dk, truth_key, by = "site", all.x = TRUE)
    dk[is.na(class), class := "not_in_truth"]
    boot <- dk[, {
      sites <- site
      cls   <- class
      n <- length(sites)
      if (n < 5L) {
        data.table(precision_se = NA_real_, recall_se = NA_real_,
                   f1_se = NA_real_, n_reps = 0L)
      } else {
        pb <- pr <- pf <- numeric(B)
        for (b in seq_len(B)) {
          ii <- sample.int(n, n, replace = TRUE)
          tr <- sum(cls[ii] == "TRUE")
          fp <- sum(cls[ii] %in% fp_classes)
          ev <- tr + fp
          pb[b] <- if (ev > 0) tr / ev else NA_real_
          pr[b] <- tr / n_true_total
          pf[b] <- if (is.finite(pb[b]) && (pb[b] + pr[b]) > 0)
                     2 * pb[b] * pr[b] / (pb[b] + pr[b]) else NA_real_
        }
        data.table(
          precision_se = sd(pb, na.rm = TRUE),
          recall_se    = sd(pr, na.rm = TRUE),
          f1_se        = sd(pf, na.rm = TRUE),
          n_reps       = as.integer(B))
      }
    }, by = gdims]
    boot[, `:=`(tier = tier_name, se_method = "site_bootstrap",
                precision_mean = NA_real_, recall_mean = NA_real_,
                f1_mean = NA_real_)]
    se_rows[[length(se_rows) + 1L]] <- boot
  } else if (tier_name == "celltype" && "arm" %in% names(d)) {
    # Cell-type tier: each cell type is an independent biological unit. Extract
    # the celltype from the arm column (e.g. AAVS1_ctrl__Adipogenic_MSC_Fib →
    # Adipogenic_MSC_Fib), score each celltype separately against the fixed
    # truth set, then report mean ± SD across cell types.
    d[, celltype := sub("^.*__", "", arm)]
    cts <- sort(unique(d$celltype))
    msg("  ", length(cts), " cell types; scoring each independently")
    ct_reps <- list()
    for (ct in cts) {
      msg("  celltype ", ct)
      dx <- unique(d[celltype == ct, c(gdims, "site"), with = FALSE])
      sc <- dx[, {
        r <- score_sites(site, truth_key, n_true_total)
        as.data.table(r)
      }, by = gdims]
      sc[, `:=`(tier = tier_name, celltype = ct)]
      ct_reps[[length(ct_reps) + 1L]] <- sc
    }
    ct_scores <- rbindlist(ct_reps, fill = TRUE)
    rep_rows[[length(rep_rows) + 1L]] <- ct_scores
    sum_sd <- ct_scores[, {
      .(n_reps         = .N,
        precision_mean = mean(precision, na.rm = TRUE),
        precision_se   = sd(precision,   na.rm = TRUE),
        recall_mean    = mean(recall_vs_union, na.rm = TRUE),
        recall_se      = sd(recall_vs_union,   na.rm = TRUE),
        f1_mean        = mean(f1, na.rm = TRUE),
        f1_se          = sd(f1,   na.rm = TRUE))
    }, by = gdims]
    sum_sd[, `:=`(tier = tier_name, se_method = "celltype_sd")]
    se_rows[[length(se_rows) + 1L]] <- sum_sd

  } else {
    msg("  no SE method for tier ", tier_name, "; skipping")
  }
}

if (length(rep_rows)) {
  reps_all <- rbindlist(rep_rows, fill = TRUE)
  # Drop accidental double-count of sample if loop appended twice
  reps_all <- unique(reps_all)
  fwrite(reps_all, file.path(outdir, "gt_tier_replicate_scores.tsv"), sep = "\t")
  msg("wrote gt_tier_replicate_scores.tsv (", nrow(reps_all), " rows)")
}

se_all <- rbindlist(se_rows, fill = TRUE)
keyc <- intersect(c("tier", "caller", "bam_type", "strand_mode", "snp_subtract"),
                  names(pooled))
out <- merge(pooled, se_all, by = keyc, all.x = TRUE)
# For sample tier, prefer replicate mean as the plotted centre when present;
# keep pooled columns as precision/recall_vs_union/f1 for continuity with
# existing figures, and expose *_mean / *_se alongside.
fwrite(out, file.path(outdir, "gt_tier_scoring_with_se.tsv"), sep = "\t")
msg("wrote gt_tier_scoring_with_se.tsv (", nrow(out), " rows)")
msg("done")
