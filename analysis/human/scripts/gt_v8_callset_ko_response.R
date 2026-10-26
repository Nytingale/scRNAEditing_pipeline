#!/usr/bin/env Rscript
# ============================================================================
# gt_v8_callset_ko_response.R   (2026-10-18)
# KO response of whole call sets (not only truth sites), merged NTC arm tier, 58 ranked configurations,
# call states none (step 4) and mismatch (step 5, main).
# For every called site (all 12 RNA mismatch types, caller strand after stranding) the stranded merged NTC and
# ADAR1 / ADAR2 KO pileups (v6 bq30, UMI-deduplicated, the ground-truth pileups) give coverage and variant reads
# (counts_for() of gt_spectrum12_eval_v6.R). A call is evaluable in a KO if m = KO cov x NTC variant fraction >= 3
# and NTC cov >= 10; it responds if its KO variant fraction <= 0.5 x NTC (t50) in an evaluable KO (either KO).
# Per call set: calls, evaluable calls, responding calls, by A-to-G vs other types.
#   resp_share_all = responding / evaluable over all calls (non-A-to-G calls respond only by chance)
#   resp_share_ag  = responding / evaluable over A-to-G calls
# Decoy chance rate = response share of non-A-to-G, non-T-to-C calls pooled over call sets.
# Positions in dbSNP common are absent from the pileup index (positions_noSNP.bed) and count as not evaluable.
# Output: human/pm/gt_eval_v8/callset_ko/callset_ko_response.tsv
# ============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("SCRNA_DIR", ".")
SP <- file.path(ROOT, "human/pm/gt_eval_v6/spectrum12")
RES <- file.path(ROOT, "human/pm/results_gt_v8/NTC_ctrl/04_filtered")
V8 <- file.path(ROOT, "human/pm/gt_eval_v8")
OUT <- file.path(V8, "callset_ko"); dir.create(OUT, showWarnings = FALSE)
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
STATE_RE <- c(none = "_step4_cov_edit_filtered\\.tsv\\.gz$", mismatch = "_step5_nomm\\.tsv\\.gz$")
M_EXP <- 3; MIN_COV <- 10; T_DROP <- 0.5
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

ranked <- unique(fread(file.path(V8, "state_summary/config_ranks.tsv"))[, c(COMBO, "config"), with = FALSE])
stopifnot(nrow(ranked) == 58L)

# ---- call sets ----------------------------------------------------------------
fs <- list.files(RES, recursive = TRUE, full.names = TRUE, pattern = "_step(4_cov_edit_filtered|5_nomm)\\.tsv\\.gz$")
meta <- data.table(f = fs)
meta[, rel := sub(paste0(RES, "/"), "", f, fixed = TRUE)]
meta[, c("caller", "bam_type", "strand_mode", "snp_dir") := tstrsplit(rel, "/", fixed = TRUE, keep = 1:4)]
meta[, snp_subtract := sub("snp_", "", snp_dir)]
meta[, state := fifelse(grepl(STATE_RE[["none"]], f), "none", "mismatch")]
meta[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
meta <- meta[config %in% ranked$config]
stopifnot(meta[, .N, by = state][, all(N == 58L)])
calls <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  x <- fread(meta$f[i], select = c("chr", "pos", "varType", "strand"), showProgress = FALSE)
  x[, `:=`(config = meta$config[i], state = meta$state[i])]
}))
msg("calls:", nrow(calls), "; unique (pos, strand, type):", uniqueN(calls[, .(chr, pos, strand, varType)]))

# ---- pileups at called positions ------------------------------------------------
pos <- fread(file.path(SP, "positions_noSNP.bed"), header = FALSE, col.names = c("chr", "s", "pos"))[, s := NULL]
need <- unique(calls[, .(chr, pos)])[pos, on = .(chr, pos), nomatch = NULL]
need[, idx := .I]
msg("called positions in the pileup index:", nrow(need), "of", uniqueN(calls[, .(chr, pos)]))
cols <- c("ref_plus", "ref_minus", "A", "C", "G", "T", "a", "c", "g", "t")
read_unit <- function(u) {
  f <- list.files(file.path(SP, "pileups"), pattern = paste0("^", u, "__.*\\.tsv\\.gz$"), full.names = TRUE)
  d <- rbindlist(lapply(f, function(x) fread(cmd = paste("gzip -cd", shQuote(x)), sep = "\t", showProgress = FALSE)))
  d <- need[d, on = .(chr, pos), nomatch = NULL]
  m <- matrix(0L, nrow(need), length(cols), dimnames = list(NULL, cols))
  for (cn in cols) m[d$idx, cn] <- as.integer(d[[cn]])
  msg(u, ":", nrow(d), "positions")
  m
}
UNITS <- c(ntc = "merged_NTC_ctrl", ko1 = "merged_ADAR1_KO", ko2 = "merged_ADARB1_KO")
U <- lapply(UNITS, read_unit)
COMP <- c(A = "T", C = "G", G = "C", T = "A")

sites <- unique(calls[, .(chr, pos, strand, varType)])
sites <- need[sites, on = .(chr, pos)]
sites[, `:=`(Y = substr(varType, 2, 2), plus = strand == "+")]
ok <- !is.na(sites$idx)
counts <- function(m) {
  alt <- refn <- rep(NA_integer_, nrow(sites))
  i <- sites$idx[ok]; pl <- sites$plus[ok]; Y <- sites$Y[ok]
  acol <- ifelse(pl, match(Y, cols), match(tolower(COMP[Y]), cols))
  alt[ok] <- m[cbind(i, acol)]
  refn[ok] <- ifelse(pl, m[i, "ref_plus"], m[i, "ref_minus"])
  list(alt = alt, cov = refn + alt)
}
S <- lapply(U, counts)
fr <- function(s) fifelse(s$cov > 0, s$alt / pmax(s$cov, 1), NA_real_)
sites[, `:=`(ntc_cov = S$ntc$cov, ntc_alt = S$ntc$alt, ntc_vf = fr(S$ntc),
             ko1_cov = S$ko1$cov, ko1_vf = fr(S$ko1), ko2_cov = S$ko2$cov, ko2_vf = fr(S$ko2))]
ev <- function(kc) !is.na(kc) & sites$ntc_cov >= MIN_COV & sites$ntc_alt > 0 & kc * sites$ntc_vf >= M_EXP
sites[, `:=`(ev1 = ev(ko1_cov), ev2 = ev(ko2_cov))]
sites[, `:=`(evaluable = ev1 | ev2,
             responds = (ev1 & ko1_vf <= (1 - T_DROP) * ntc_vf) | (ev2 & ko2_vf <= (1 - T_DROP) * ntc_vf))]
sites[is.na(responds), responds := FALSE]
for (t in c(0.25, 0.75)) {
  rc <- paste0("responds_t", 100 * t)
  sites[, (rc) := (ev1 & ko1_vf <= (1 - t) * ntc_vf) | (ev2 & ko2_vf <= (1 - t) * ntc_vf)]
  sites[is.na(get(rc)), (rc) := FALSE]
}
sites[, class := fifelse(varType == "AG", "AG", fifelse(varType == "TC", "TC", "decoy"))]
msg("evaluable sites by class:"); print(sites[, .(n = .N, evaluable = sum(evaluable), responds = sum(responds)), by = class])

calls <- sites[, .(chr, pos, strand, varType, class, evaluable, responds, responds_t25, responds_t75)][
  calls, on = .(chr, pos, strand, varType)]

# thresholds t25 / t50 / t75 (as the KO truth sets)
thr <- calls[, .(n_eval = sum(evaluable), n_resp_t25 = sum(responds_t25), n_resp_t50 = sum(responds),
                 n_resp_t75 = sum(responds_t75),
                 n_eval_decoy = sum(evaluable & class == "decoy"),
                 n_resp_decoy_t25 = sum(responds_t25 & class == "decoy"),
                 n_resp_decoy_t50 = sum(responds & class == "decoy"),
                 n_resp_decoy_t75 = sum(responds_t75 & class == "decoy"),
                 n_eval_ag = sum(evaluable & class == "AG"),
                 n_resp_ag_t25 = sum(responds_t25 & class == "AG"), n_resp_ag_t50 = sum(responds & class == "AG"),
                 n_resp_ag_t75 = sum(responds_t75 & class == "AG")), by = .(state, config)]
fwrite(thr, file.path(OUT, "callset_ko_response_thresholds.tsv"), sep = "\t")

# calls removed by the step-5 filter (in step 4, absent from step 5) vs calls kept
key <- c("config", "chr", "pos", "strand", "varType")
kept <- unique(calls[state == "mismatch", ..key])[, kept := TRUE]
rm4 <- kept[calls[state == "none"], on = key][is.na(kept), kept := FALSE]
rem <- rm4[, .(n_calls = .N, n_ag = sum(class == "AG"), n_eval = sum(evaluable),
               n_resp_t25 = sum(responds_t25), n_resp_t50 = sum(responds), n_resp_t75 = sum(responds_t75),
               n_eval_ag = sum(evaluable & class == "AG"), n_resp_ag_t25 = sum(responds_t25 & class == "AG"),
               n_resp_ag_t50 = sum(responds & class == "AG"), n_resp_ag_t75 = sum(responds_t75 & class == "AG")),
           by = .(config, kept)]
fwrite(rem, file.path(OUT, "callset_ko_removed_vs_kept.tsv"), sep = "\t")
res <- calls[, .(n_calls = .N, n_ag = sum(class == "AG"),
                 n_eval = sum(evaluable), n_resp = sum(responds),
                 n_eval_ag = sum(evaluable & class == "AG"), n_resp_ag = sum(responds & class == "AG"),
                 n_eval_decoy = sum(evaluable & class == "decoy"), n_resp_decoy = sum(responds & class == "decoy")),
             by = .(state, config)]
res[, `:=`(resp_share_all = n_resp / n_eval, resp_share_ag = n_resp_ag / n_eval_ag)]
chance <- res[, sum(n_resp_decoy) / sum(n_eval_decoy)]
res[, decoy_chance := chance]
res <- merge(ranked, res, by = "config")
fwrite(res, file.path(OUT, "callset_ko_response.tsv"), sep = "\t")
msg("decoy chance response share:", round(chance, 4))
print(res[, .(med_all = median(resp_share_all), med_ag = median(resp_share_ag)), by = state])
