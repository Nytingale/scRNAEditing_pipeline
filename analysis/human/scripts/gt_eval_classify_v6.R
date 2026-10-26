#!/usr/bin/env Rscript
# =============================================================================
# gt_eval_classify_v6.R   (copy of gt_eval_classify_v3.R used for GT v6)
#
# v6 = genome-wide, caller-independent truth set (human/pm/gt_eval_v6/,
# driver human/pm/scripts/run_gt_truth_rebuild_v6.sh). Criteria logic is
# byte-for-byte v3. The ONLY change: --rediportal_sites <file of chr:pos>
# supplies REDIportal membership pre-computed by bedtools intersect, replacing
# the per-run `bcftools query -R` lookup (too slow on a genome-wide BED).
# When --rediportal_sites is absent the v3 bcftools path is used unchanged.
# 2026-09-27: --min_cov_teratoma (per-teratoma evaluability for consistency)
# and --min_cov_ko (KO-arm depth), both default = --min_cov, i.e. identical to
# v3 unless set.
# 2026-09-27: --ko_min_expected m (KO arm counts only if KO cov x control
# editing >= m, a depth-invariant informative-KO gate); default 0 = off.
#
# ----------------------------- original v3 header -----------------------------
# gt_eval_classify_v3.R   (dataset-agnostic GT evaluator, step 3/3 -- v3)
#
# Implements the v3 ground-truth definition
# (human/pm/gt_eval/true_def_v3/GROUND_TRUTH_DEFINITION.md) on a TIER-NEUTRAL
# architecture.
#
# ---------------------------------------------------------------------------
# WHY THIS REPLACES gt_eval_classify.R RATHER THAN EXTENDING IT
#
# The v1 engine (still on disk, untouched) builds its reference set FROM the
# combination being scored:  "sites this combo discovered, covered in both
# arms". That couples the yardstick to the thing being measured, and it has
# three consequences the v3 definition cannot live with:
#   * truth is re-derived per tier, so a cell-type tier must establish its own
#     reproducibility -- impossible, because per-celltype BAMs pool T1-T4
#     (HANDOFF_2026-09-16_pm §4c objection 1);
#   * a site's label can differ between tiers, so tier precisions are not
#     comparable and §4c had to recommend never printing them side by side;
#   * it silently inherits whatever site list the pileups were built at -- the
#     bug found 2026-09-17, where pileups were at the MERGED-arm union so
#     35,241 per-sample-discovered sites (2,717 of them scallele-exclusive)
#     were invisible.
#
# v3 instead computes truth ONCE, over the union of every tier's discoveries,
# from pileups that touch neither the callers nor the tiers. A site gets ONE
# label. Tiers are then scored against that fixed yardstick, which makes them
# directly comparable and dissolves all three problems above.
#
# v5 (--site_list per tier) builds a SEPARATE truth set for each tier: Phase 1
# classifies only that tier's discovery sites. No cross-tier union enters the
# truth table or the recall denominator. Pileups are still shared; the site
# universe is not.
#
# ---------------------------------------------------------------------------
# THE DEFINITION (§2a), applied in this order
#   1. >= --min_cov reads at the site in control AND in the KO arm
#   2. control edited: >= --min_edit_reads edited reads AND >= --min_ctrl_edit
#   3. control-positive in >= --consistency_min of the teratomas
#   4. depleted: KO <= (1 - d) * control.  Primary d = --depletion
#   5. depleted in ADAR1_KO *OR* ADARB1_KO.  ADARB2/ADAR3 excluded entirely
#   6. four classes: TRUE / FP_novel / FP_known_not_depleted / ambiguous
#   7. no differential-editing test -- effect size + reproducibility only
#
# Magnitude (1,2,4) is measured on the MERGED arms, the deepest available view;
# reproducibility (3) on the per-teratoma BAMs, the only place replicate
# structure exists. Both read the SAME combined site list, which is what makes
# the label tier-independent.
#
# ---------------------------------------------------------------------------
# TWO DEPARTURES FROM A LITERAL READING OF §2a, BOTH DELIBERATE
#
# (a) --consistency_mode evaluable  (default; `strict` reproduces §2a exactly)
#     Per-unit depth varies ~3-5x WITHIN an arm -- measured 2026-09-17 from
#     human/pm/gt_eval/tier_depth.tsv:
#         T4_ADAR1_KO  1,126,629 reads   vs  T2_ADAR1_KO  5,480,579  (4.9x)
#         T2_NTC_ctrl 19,895,031 reads   vs  T4_NTC_ctrl 66,024,598  (3.3x)
#     A literal "3 of 4" therefore fails sites purely because one teratoma is
#     shallow there, which is a depth artefact reported as biological
#     non-reproducibility. `evaluable` mode requires 3 of 4 among the teratomas
#     that actually reach --min_cov, and records n_evaluable so the
#     distinction stays visible. Both modes are written out.
#
# (b) REDIportal is NOT folded into precision (§2f). It splits the FP class
#     into FP_known_not_depleted vs FP_novel and is reported on its own axis.
#     85.7% of the discovery set is already catalogued, so using it as an
#     adjustment absolves FPs almost mechanically.
#
# ---------------------------------------------------------------------------
# OUTPUTS (in --outdir)
#   gt_truth_sites.tsv.gz    one row per site: class, powered, per-arm levels
#   gt_truth_summary.tsv     class counts across the depletion series
#   gt_tier_scoring.tsv      per tier x caller x bam_type x strand_mode
#   (plots are produced by the BENCHMARK_* scripts, not here)
#
# USAGE
#   # phase 1 only -- build the truth set
#   Rscript gt_eval_classify_v3.R --truth_pileup_dir human/pm/gt_eval/truth_combined \
#       --outdir human/pm/gt_eval/true_def_v3
#
#   # phase 1 + 2 -- also score whichever tiers exist
#   Rscript gt_eval_classify_v3.R --truth_pileup_dir ... --outdir ... \
#       --discovery "sample=.../persample_union/gt_eval_discovery.tsv.gz,dataset=.../v2/gt_eval_discovery.tsv.gz"
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x)) || !length(x)) y else x
parse_cli <- function() {
  a <- commandArgs(trailingOnly = TRUE); out <- list(); i <- 1L
  while (i <= length(a)) {
    if (grepl("^--", a[i])) { out[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L } else i <- i + 1L
  }
  out
}
args <- parse_cli()

pile_dir  <- args[["truth_pileup_dir"]] %||% stop("--truth_pileup_dir required")
outdir    <- args[["outdir"]]           %||% stop("--outdir required")
ctrl_arm  <- args[["control_arm"]]      %||% "NTC_ctrl"
treat_arms<- strsplit(args[["treatment_arms"]] %||% "ADAR1_KO,ADARB1_KO", ",")[[1]]
teratomas <- strsplit(args[["teratomas"]]      %||% "T1,T2,T3,T4", ",")[[1]]
min_cov        <- as.numeric(args[["min_cov"]]        %||% 10)
# v6: per-teratoma evaluability depth for the consistency rule, decoupled from
# the merged-arm --min_cov (which also gates the KO arms). Defaults to
# --min_cov, so omitting it reproduces v3 exactly.
min_cov_teratoma <- as.numeric(args[["min_cov_teratoma"]] %||% min_cov)
# v6: KO-arm depth, decoupled from the merged-control --min_cov (the KO arms are
# far shallower). Defaults to --min_cov, so omitting it reproduces v3 exactly.
min_cov_ko <- as.numeric(args[["min_cov_ko"]] %||% min_cov)
# v6: depth-invariant informative-KO gate. A KO arm counts only if the edited
# reads it would carry under NO depletion (KO cov x merged-control editing) are
# at least this many. 0 (default) disables it, reproducing earlier runs exactly.
ko_min_expected <- as.numeric(args[["ko_min_expected"]] %||% 0)
min_edit_reads <- as.numeric(args[["min_edit_reads"]] %||% 3)
min_ctrl_edit  <- as.numeric(args[["min_ctrl_edit"]]  %||% 0.05)
consistency_min  <- as.numeric(args[["consistency_min"]]  %||% 3)
consistency_mode <- args[["consistency_mode"]] %||% "evaluable"   # evaluable | strict
depletion      <- as.numeric(args[["depletion"]] %||% 0.50)
# Default series is a 10% grid UNION the levels §2a.4 commits to reporting
# (25% and 75%), so the curve is evenly sampled without dropping the points the
# definition names. Every level is computed from scratch by class_at(), so
# nothing here is interpolated.
dep_series     <- as.numeric(strsplit(
  args[["depletion_series"]] %||%
    "0.10,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.75,0.80,0.90", ",")[[1]])
ci_alpha       <- as.numeric(args[["ci_alpha"]] %||% 0.05)
rediportal     <- args[["rediportal"]] %||% NA_character_
union_bed      <- args[["union_bed"]]  %||% NA_character_
redi_sites_arg <- args[["rediportal_sites"]] %||% NA_character_
disc_spec      <- args[["discovery"]]  %||% NA_character_
dataset_label  <- args[["dataset_label"]] %||% "Teratoma_KO"
# Optional site universe restriction. When set, PHASE 1 classifies ONLY these
# sites (the tier's own discovery list for v5). This is a separate truth set,
# not a post-hoc filter on a cross-tier union classification.
site_list_arg  <- args[["site_list"]] %||% NA_character_
# truth_scope retained for backwards compatibility; v5 uses --site_list per tier.
truth_scope    <- args[["truth_scope"]] %||% "union"
if (!truth_scope %in% c("union", "per_tier"))
  stop("--truth_scope must be 'union' or 'per_tier', got: ", truth_scope)

# `none` drops criterion 3 entirely, for datasets with NO replicate structure.
# Replogle is one K562 cell line with three arms and no donors, so each arm is a
# single giant sample: there is nothing to be reproducible ACROSS, and a
# reproducibility filter would either be vacuous or fabricate a replicate axis
# out of the arms (which are the contrast, not replicates of each other).
#
# >>> THIS MAKES THE TRUTH SETS DIFFERENTLY CONSTRUCTED. <<<
# Teratoma's truth set IS reproducibility-filtered and Replogle's cannot be, so
# a cross-dataset comparison is not like-for-like on that axis. The comparable
# teratoma arm is this same `none` mode, measured 2026-09-18 at d = 0.50:
#
#   consistency rule            TRUE (point)   powered
#   evaluable (adopted)                4,754     1,063
#   none (Replogle-comparable)        21,560     1,165
#
# Note `none` is MORE permissive than the "n_positive >= 1" arm in the
# consistency sensitivity table (17,181 / 1,158): that arm still demands the
# site be control-positive in at least one teratoma, i.e. per-teratoma coverage
# AND edited reads AND >=5%, whereas `none` asks only that the MERGED control be
# edited. They are different relaxations and must not be quoted for each other.
#
# The headline is unaffected either way: the powered core moves 1,063 -> 1,165
# (10%) while the nominal count moves 4.5-fold. Report teratoma in `none` mode
# alongside Replogle, and the adopted rule for teratoma on its own.
stopifnot(consistency_mode %in% c("evaluable", "strict", "none"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

# ---------------------------------------------------------------------------
# Pileup reader. parse_pileup_stranded.pl emits, per covered position:
#   chr pos ref coverage variant_cov A C G T a c g t
#   ref_count_plus ref_count_minus editing_level
# It is already strand-aware (direction derived from the genomic ref base and
# the library type), so no strand handling is repeated here. variant_cov IS the
# edited-read count, and editing_level == variant_cov / coverage -- verified
# 2026-09-17. Positions with zero coverage are simply absent from the file and
# are treated as cov 0.
read_pileup <- function(unit) {
  f <- file.path(pile_dir, paste0("pileup_", unit, ".tsv"))
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  p <- fread(f, sep = "\t", header = TRUE, showProgress = FALSE)
  if (!nrow(p)) return(NULL)
  p[, site := paste0(chr, ":", pos)]
  p[, .(site,
        cov   = as.numeric(coverage),
        alt   = as.numeric(variant_cov),
        edit  = as.numeric(editing_level))]
}

# ===========================================================================
# PHASE 1 -- BUILD THE TRUTH SET
# ===========================================================================
msg("PHASE 1: building tier-neutral truth set")
msg("  control=", ctrl_arm, "  treatments=", paste(treat_arms, collapse = ","),
    "  (ADARB2/ADAR3 excluded by definition §2a.5)")

# --- magnitude, from the merged arms ---------------------------------------
merged <- list()
for (arm in c(ctrl_arm, treat_arms)) {
  p <- read_pileup(paste0("merged_", arm))
  if (is.null(p)) stop("missing or empty merged pileup for arm '", arm,
                       "' in ", pile_dir)
  setnames(p, c("cov", "alt", "edit"),
           paste0(c("cov_", "alt_", "edit_"), arm))
  merged[[arm]] <- p
  msg("  merged_", arm, ": ", nrow(p), " covered sites")
}
truth <- Reduce(function(a, b) merge(a, b, by = "site", all = TRUE), merged)
for (j in setdiff(names(truth), "site"))
  set(truth, which(is.na(truth[[j]])), j, 0)
msg("  sites in merged pileup union: ", nrow(truth))

# --- optional site-universe restriction (v5 per-tier truth) -----------------
# Applied BEFORE criteria 1–5 so the truth set is built only over this tier's
# discoveries. Labels for overlapping sites match a union run (criteria are
# site-local), but n_true_total / gt_truth_summary belong to this universe
# alone -- no cross-tier sites enter the denominator or the truth table.
if (!is.na(site_list_arg) && nzchar(site_list_arg)) {
  if (!file.exists(site_list_arg))
    stop("--site_list not found: ", site_list_arg)
  keep <- unique(fread(site_list_arg, header = FALSE, col.names = "site",
                       showProgress = FALSE)$site)
  n_before <- nrow(truth)
  truth <- truth[site %in% keep]
  msg("  --site_list: kept ", nrow(truth), " / ", n_before,
      " pileup sites (list has ", length(keep), " sites; ",
      length(keep) - nrow(truth), " discovery sites have no pileup row)")
}

# --- reproducibility, from the per-teratoma control BAMs --------------------
# Criterion 3 is about the CONTROL arm: is the site reliably edited in control
# across teratomas. KO arms are not used here.
ctrl_pos <- NULL; ctrl_eval <- NULL
if (consistency_mode == "none") {
  msg("  consistency DROPPED (--consistency_mode none): no replicate structure")
  msg("  every arm is treated as one sample; criteria 1, 2, 4 and 5 still apply")
} else
for (tt in teratomas) {
  p <- read_pileup(paste0(tt, "_", ctrl_arm))
  if (is.null(p)) { msg("  WARNING: no pileup for ", tt, "_", ctrl_arm, "; skipped"); next }
  p[, evaluable := cov >= min_cov_teratoma]
  p[, positive  := cov >= min_cov_teratoma & alt >= min_edit_reads & edit >= min_ctrl_edit]
  e <- p[evaluable == TRUE, .(site)][, (tt) := 1L]
  o <- p[positive  == TRUE, .(site)][, (tt) := 1L]
  ctrl_eval <- if (is.null(ctrl_eval)) e else merge(ctrl_eval, e, by = "site", all = TRUE)
  ctrl_pos  <- if (is.null(ctrl_pos))  o else merge(ctrl_pos,  o, by = "site", all = TRUE)
  msg("  ", tt, "_", ctrl_arm, ": ", nrow(e), " evaluable, ", nrow(o), " control-positive")
}
if (consistency_mode == "none") {
  # Columns are still written so the output schema does not change between
  # datasets; they are explicitly NA rather than 0, because "not assessable"
  # and "assessed and found in zero samples" are different statements.
  truth[, `:=`(n_positive = NA_integer_, n_evaluable = NA_integer_,
               consistent_strict = NA, consistent_evaluable = NA,
               consistent = TRUE)]
} else {
tcols <- intersect(teratomas, names(ctrl_pos))
for (j in tcols) {
  set(ctrl_pos,  which(is.na(ctrl_pos[[j]])),  j, 0L)
  set(ctrl_eval, which(is.na(ctrl_eval[[j]])), j, 0L)
}
ctrl_pos[,  n_positive  := rowSums(.SD), .SDcols = tcols]
ctrl_eval[, n_evaluable := rowSums(.SD), .SDcols = tcols]
# Per-teratoma flags are carried through, not just their totals: the
# depth-confounding panel needs positivity PER teratoma against that
# teratoma's own depth, which a count cannot show.
setnames(ctrl_pos,  tcols, paste0("pos_",  tcols))
setnames(ctrl_eval, tcols, paste0("eval_", tcols))
truth <- merge(truth, ctrl_pos,  by = "site", all.x = TRUE)
truth <- merge(truth, ctrl_eval, by = "site", all.x = TRUE)
for (j in c("n_positive", "n_evaluable",
            paste0("pos_", tcols), paste0("eval_", tcols)))
  set(truth, which(is.na(truth[[j]])), j, 0L)

# strict   : >= k of ALL teratomas (literal §2a.3) -- fails a site when one
#            teratoma merely lacks depth there.
# evaluable: "positive in AT LEAST TWO independent teratomas, with NO MORE THAN
#            ONE DROPOUT among the teratomas where the site is evaluable", i.e.
#              need = max(2, n_evaluable - 1)
#              2 evaluable -> 2 needed     3 evaluable -> 2 needed
#              4 evaluable -> 3 needed     (identical to strict at n = 4)
#            Two clauses, each doing a distinct job: the `2` floor means one
#            deep teratoma can never establish reproducibility on its own; the
#            `n_evaluable - 1` clause means a site is never failed for lacking
#            depth in a teratoma rather than lacking signal.
#
#            This REPLACES an earlier proportional form, pmax(2, round(0.75 *
#            n_evaluable)). The two are arithmetically IDENTICAL for
#            n_evaluable <= 5 and so change nothing on teratoma's four, but
#            they diverge from n_evaluable >= 6 (proportional asks 4 of 6 where
#            this asks 5 of 6). The explicit form is kept because it states the
#            criterion the way it is reported, and because it generalises
#            correctly to datasets with more samples.
truth[, consistent_strict := n_positive >= consistency_min]
truth[, consistent_evaluable :=
        n_positive >= 2 &
        n_positive >= pmax(2, n_evaluable - 1)]
truth[, consistent := if (consistency_mode == "strict") consistent_strict
                      else consistent_evaluable]
}

# --- criteria 1 and 2, on the merged arms -----------------------------------
ctrl_cov  <- truth[[paste0("cov_",  ctrl_arm)]]
ctrl_alt  <- truth[[paste0("alt_",  ctrl_arm)]]
ctrl_edit <- truth[[paste0("edit_", ctrl_arm)]]
truth[, ctrl_edited := ctrl_alt >= min_edit_reads & ctrl_edit >= min_ctrl_edit]

# --- criteria 4 and 5, per depletion level ----------------------------------
# Clopper-Pearson one-sided upper limit on the KO rate: qbeta(1-a, k+1, n-k).
# `powered` is a DEPTH property: even with zero edited reads observed, is n big
# enough for the bound to fall below the target (at k=0 the bound is
# 1 - a^(1/n))?  `true_ci` is the same bound evaluated at the OBSERVED k.
# When every read is edited (k == n) the upper limit is 1 by definition, and
# qbeta's second shape parameter would be 0; those rows keep the initialised 1.
cp_upper <- function(k, n, a = ci_alpha) {
  out <- rep(1, length(k))
  ok <- n > 0 & (n - k) > 0
  out[ok] <- qbeta(1 - a, k[ok] + 1, n[ok] - k[ok])
  out
}

class_at <- function(d) {
  target <- (1 - d) * ctrl_edit
  dep_any <- rep(FALSE, nrow(truth))
  pow_any <- rep(FALSE, nrow(truth))
  ci_any  <- rep(FALSE, nrow(truth))
  cov_any <- rep(FALSE, nrow(truth))
  for (arm in treat_arms) {
    kc <- truth[[paste0("cov_",  arm)]]
    ka <- truth[[paste0("alt_",  arm)]]
    ke <- truth[[paste0("edit_", arm)]]
    covered <- ctrl_cov >= min_cov & kc >= min_cov_ko & kc * ctrl_edit >= ko_min_expected
    cov_any <- cov_any | covered
    dep_any <- dep_any | (covered & ke <= target)
    pow_any <- pow_any | (covered & cp_upper(rep(0, length(kc)), kc) < target)
    ci_any  <- ci_any  | (covered & cp_upper(ka, kc) < target)
  }
  list(covered = cov_any, depleted = dep_any, powered = pow_any, ci = ci_any)
}

# --- REDIportal membership (own axis, §2f) ----------------------------------
in_redi <- rep(NA, nrow(truth))
if (!is.na(redi_sites_arg) && nzchar(redi_sites_arg)) {
  if (!file.exists(redi_sites_arg))
    stop("--rediportal_sites not found: ", redi_sites_arg)
  redi_sites <- unique(fread(redi_sites_arg, header = FALSE, col.names = "site",
                             showProgress = FALSE)$site)
  in_redi <- truth$site %in% redi_sites
  msg("  REDIportal (pre-computed, ", length(redi_sites), " listed): ",
      sum(in_redi), " of ", nrow(truth), " sites (",
      round(100 * mean(in_redi), 1), "%)")
} else if (!is.na(rediportal) && file.exists(rediportal) && !is.na(union_bed) && file.exists(union_bed)) {
  msg("  querying REDIportal membership")
  cmd <- paste("bcftools query -f '%CHROM\\t%POS\\n' -R", shQuote(union_bed),
               shQuote(rediportal), "2>/dev/null")
  r <- tryCatch(fread(cmd = cmd, header = FALSE, col.names = c("chr", "pos"),
                      showProgress = FALSE),
                error = function(e) data.table())
  if (nrow(r)) {
    redi_sites <- unique(paste0(r$chr, ":", r$pos))
    in_redi <- truth$site %in% redi_sites
    msg("  REDIportal: ", sum(in_redi), " of ", nrow(truth), " sites (",
        round(100 * mean(in_redi), 1), "%)")
  } else msg("  WARNING: REDIportal query returned nothing; membership left NA")
} else {
  msg("  REDIportal not supplied; FP class will not be split (§2f axis unavailable)")
}
truth[, in_rediportal := in_redi]

# --- assemble the primary classification ------------------------------------
prim <- class_at(depletion)
truth[, evaluable := prim$covered]
truth[, powered   := prim$powered]
truth[, depleted  := prim$depleted]
truth[, depleted_ci := prim$ci]

have_redi <- any(!is.na(in_redi))
truth[, class := "ambiguous"]
# TRUE: covered, control genuinely edited, reproducible, and depleted in either arm
truth[evaluable & ctrl_edited & consistent & depleted, class := "TRUE"]
# FP: covered, control genuinely edited and reproducible, but NOT depleted anywhere.
# Split by REDIportal only when membership is actually available -- otherwise a
# missing catalogue would silently relabel every FP as novel, which is the
# opposite of the truth.
if (have_redi) {
  truth[evaluable & ctrl_edited & consistent & !depleted & in_rediportal,
        class := "FP_known_not_depleted"]
  truth[evaluable & ctrl_edited & consistent & !depleted & !in_rediportal,
        class := "FP_novel"]
} else {
  truth[evaluable & ctrl_edited & consistent & !depleted, class := "FP_unsplit"]
}
# everything else stays ambiguous: not covered, control not edited, or not reproducible
fp_classes <- c("FP_novel", "FP_known_not_depleted", "FP_unsplit")

msg("  primary classification at d=", depletion, ":")
print(truth[, .N, by = class][order(-N)])

# --- depletion series (§2c) --------------------------------------------------
ser <- rbindlist(lapply(dep_series, function(d) {
  cl <- class_at(d)
  base <- truth$evaluable & truth$ctrl_edited & truth$consistent
  n_true <- sum(base & cl$depleted)
  n_fp   <- sum(base & !cl$depleted)
  # THREE distinct "powered" quantities, kept separate on purpose because
  # conflating them is easy and misleading:
  #   n_depth_sufficient  gate AND powered            -- could the site ever show
  #                                                      depletion, whatever the
  #                                                      outcome? (a denominator)
  #   n_true_powered      gate AND depleted AND powered -- THE HEADLINE: depleted,
  #                                                      with the depth for that
  #                                                      to mean something
  #   n_true_ci           gate AND ci                  -- depletion actually proven
  # `n_powered` is retained under its old name for continuity with existing
  # tables, but it is the DENOMINATOR quantity and must not be read as a count
  # of true sites; it was previously plotted beside panel (d)'s intersection
  # under the same word, giving 3,251 and 1,063 both labelled "powered".
  data.table(depletion = d,
             n_true_point       = n_true,
             n_powered          = sum(base & cl$powered),
             n_depth_sufficient = sum(base & cl$powered),
             n_true_powered     = sum(base & cl$depleted & cl$powered),
             n_true_ci          = sum(base & cl$ci),
             n_fp               = n_fp,
             precision    = if ((n_true + n_fp) > 0) n_true / (n_true + n_fp) else NA_real_)
}))
msg("  depletion series:")
print(ser)

fwrite(truth, file.path(outdir, "gt_truth_sites.tsv.gz"), sep = "\t")
fwrite(ser,   file.path(outdir, "gt_truth_summary.tsv"),  sep = "\t")
msg("  wrote gt_truth_sites.tsv.gz (", nrow(truth), " sites) and gt_truth_summary.tsv")

# ===========================================================================
# PHASE 2 -- SCORE EACH TIER AGAINST THE TRUTH SET
# ===========================================================================
# NOTE ON RECALL (truth_scope = union, v3/v4). The truth set is defined over
# the union of what the tiers discovered, so a site nobody found cannot be in
# it. Recall is therefore relative to the CROSS-TIER UNION and is labelled
# n_true_recovered / n_true_total.
#
# NOTE ON RECALL (truth_scope = per_tier, v5). Each tier is its own truth set:
# n_true_total = count of TRUE sites among that tier's discovery sites only.
# A TRUE site found only by another tier does not inflate this tier's
# denominator. Column name recall_vs_union is kept for schema compatibility
# with figure scripts; under per_tier it means recall vs the tier truth set.
#
# NOTE ON EVALUABILITY. The merged/dataset tier discovers preferentially where
# it has depth, and evaluability also needs depth, so its discoveries are
# enriched for evaluable sites. That inflates evaluability_rate but NOT
# precision, because among evaluable sites the label is independent of who
# discovered it (the merged pileup at a site is the same reads either way).
# evaluability_rate is reported as a first-class column so the advantage is
# visible rather than hidden.
if (!is.na(disc_spec)) {
  msg("PHASE 2: scoring tiers (truth_scope=", truth_scope, ")")
  truth_key <- truth[, .(site, class, evaluable, powered)]
  n_true_global <- sum(truth$class == "TRUE")

  specs <- strsplit(strsplit(disc_spec, ",")[[1]], "=")
  rows <- list()
  tier_summ <- list()
  for (s in specs) {
    if (length(s) != 2L) { msg("  skipping malformed --discovery entry"); next }
    tier_name <- s[1]; path <- s[2]
    if (!file.exists(path)) { msg("  SKIP tier '", tier_name, "': no discovery table at ", path); next }
    d <- tryCatch(fread(cmd = paste("gzip -cd", shQuote(path)), sep = "\t", showProgress = FALSE),
                  error = function(e) data.table())
    if (!nrow(d)) { msg("  SKIP tier '", tier_name, "': empty discovery table"); next }
    msg("  tier '", tier_name, "': ", nrow(d), " discovery rows, ",
        uniqueN(d$combo), " combos")
    # Group by EVERY grid dimension present. Reporting a ranking by caller alone
    # is not answerable -- "reditools3" differs by dedup, stranding and SNP
    # subtraction, and each gives different numbers - so the unit of comparison
    # is the full unique combination. snp_subtract is absent from older
    # discovery tables, hence the intersect().
    gdims <- intersect(c("caller", "bam_type", "strand_mode", "snp_subtract"), names(d))
    if (!"snp_subtract" %in% gdims)
      msg("  NOTE tier '", tier_name, "': no snp_subtract column; ",
          "combinations are pooled over that dimension")
    disc_sites <- unique(d$site)
    if (truth_scope == "per_tier") {
      # Tier truth set = classified sites that this tier discovered. Labels
      # themselves are unchanged (site-local); only the universe shrinks.
      n_true_total <- sum(truth$site %in% disc_sites & truth$class == "TRUE")
      n_fp_tier    <- sum(truth$site %in% disc_sites & truth$class %in% fp_classes)
      n_amb_tier   <- sum(truth$site %in% disc_sites & truth$class == "ambiguous")
      tier_summ[[length(tier_summ) + 1L]] <- data.table(
        tier = tier_name,
        n_discovery_sites = length(disc_sites),
        n_true_point = n_true_total,
        n_fp = n_fp_tier,
        n_ambiguous = n_amb_tier,
        n_true_global_union = n_true_global)
      msg("    per-tier truth: n_true=", n_true_total,
          " / ", length(disc_sites), " discovery sites",
          " (global union TRUE=", n_true_global, ")")
    } else {
      n_true_total <- n_true_global
    }
    dk <- unique(d[, c(gdims, "site"), with = FALSE])
    dk <- merge(dk, truth_key, by = "site", all.x = TRUE)
    dk[is.na(class), class := "not_in_truth"]
    rows[[length(rows) + 1L]] <- dk[, {
      n_disc <- .N
      n_tr   <- sum(class == "TRUE")
      n_fp   <- sum(class %in% fp_classes)
      n_eval <- n_tr + n_fp
      .(tier = tier_name,
        n_discovered      = n_disc,
        n_evaluable       = n_eval,
        evaluability_rate = if (n_disc > 0) n_eval / n_disc else NA_real_,
        n_true            = n_tr,
        n_fp              = n_fp,
        n_fp_novel        = sum(class == "FP_novel"),
        n_fp_known        = sum(class == "FP_known_not_depleted"),
        n_ambiguous       = sum(class == "ambiguous"),
        precision         = if ((n_tr + n_fp) > 0) n_tr / (n_tr + n_fp) else NA_real_,
        n_true_recovered  = n_tr,
        n_true_total      = n_true_total,
        recall_vs_union   = if (n_true_total > 0) n_tr / n_true_total else NA_real_)
    }, by = gdims]
  }
  if (length(tier_summ)) {
    fwrite(rbindlist(tier_summ), file.path(outdir, "gt_truth_summary_by_tier.tsv"),
           sep = "\t")
    msg("  wrote gt_truth_summary_by_tier.tsv")
  }
  if (length(rows)) {
    scoring <- rbindlist(rows, fill = TRUE)
    setcolorder(scoring, intersect(c("tier", "caller", "bam_type",
                                     "strand_mode", "snp_subtract"),
                                   names(scoring)))
    setorder(scoring, tier, -precision, -n_true)
    fwrite(scoring, file.path(outdir, "gt_tier_scoring.tsv"), sep = "\t")
    msg("  wrote gt_tier_scoring.tsv (", nrow(scoring), " tier x method rows)")
    print(head(scoring, 15))
  } else {
    msg("  no tiers scored -- nothing written")
  }
} else {
  msg("PHASE 2 skipped (no --discovery given)")
}

msg("done")
