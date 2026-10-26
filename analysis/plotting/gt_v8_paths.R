# ============================================================================
# gt_v8_paths.R   (2026-10-13)
# The v8 paper truth, shared by publication/teratoma_gt/scripts/TERATOMA_common.R and
# publication/gt_definition/scripts/GTDEF_common.R. Definition: human/pm/gt_eval_v8/README.md.
# Requires REPO_DIR.
# ============================================================================
V8_DIR        <- file.path(REPO_DIR, "human", "pm", "gt_eval_v8")
V8_LABEL      <- "v8"
V8_SRC_VARIANT <- "v7_nomm_pgp1"   # the gt_eval_v7 variant v8 was copied from
V8_STATE      <- "mismatch"        # call state: nearby-mismatch pass only
V8_STATES     <- c("none", "mismatch", "clustered", "both")
V8_TRUTH_SITES <- file.path(V8_DIR, "truth", "gt_v7_truth_sites.tsv.gz")
V8_TRUTH_COUNTS <- file.path(V8_DIR, "truth", "gt_v8_truth_counts.tsv")
V8_SUPPORT    <- file.path(V8_DIR, "truth", "support")
# Recall denominator (2026-10-22): "pwcov" (headline; scoring_pwcov/, GT_UNIT_RULE=gate_cov, 8 cell types with
# >= 30 truth sites) or V8_SCORING=cov (scoring/, >= 10 reads for every truth, 11 cell types with >= 20)
V8_SCORING    <- Sys.getenv("V8_SCORING", "pwcov")
stopifnot(V8_SCORING %in% c("pwcov", "cov"))
V8_SFX        <- if (V8_SCORING == "pwcov") "_pwcov" else ""
V8_CMP_RANK   <- file.path(V8_DIR, paste0("compare_config_rank", V8_SFX, ".tsv"))
v8_scoring_dir <- function(state = V8_STATE) file.path(V8_DIR, paste0("scoring", V8_SFX), state)
v8_config_rank <- function(state = V8_STATE) file.path(v8_scoring_dir(state), "config_rank.tsv")
v8_redi_dir    <- function(state = V8_STATE) file.path(V8_DIR, paste0("rediportal", V8_SFX), paste0("scoring_", state))
V8_REDI_VARIANTS <- c("cov10_alt1", "alt3", "alt3_e05")
v8_redi_files  <- function(state = V8_STATE)
  setNames(file.path(v8_redi_dir(state), paste0("gt_tier_scoring_", V8_REDI_VARIANTS, ".tsv")),
           paste0("redi_", V8_REDI_VARIANTS))
V8_CT_KEEP    <- file.path(V8_DIR, "unit_cov",                                 # cell-type tier units
                           if (V8_SCORING == "pwcov") "celltype_keep_pwcov.txt" else "celltype_keep_min20.txt")
V8_CT_COVERED <- file.path(V8_DIR, "unit_cov",
                           if (V8_SCORING == "pwcov") "celltype_covered_positives_pwcov.tsv" else "celltype_covered_positives.tsv")
V8_CT_MIN_SITES <- if (V8_SCORING == "pwcov") 30L else 20L
V8_N_CELLTYPES <- length(readLines(V8_CT_KEEP))
V8_KO_TRUTHS  <- c("expected_t25", "expected_t50", "expected_t75", "depth_t25", "depth_t50", "depth_t75")
V8_TRUTH_IDS  <- c(V8_KO_TRUTHS, paste0("redi_", V8_REDI_VARIANTS))
# Figure text: "v8" and pipeline names never appear in a figure.
V8_FILTER_LAB <- "Removed if near other mismatch types"
# Gate names: power-gated = >= 3 expected edited KO reads (frozen, headline); depth-gated = >= 10 KO reads
V8_GATE_LAB   <- c(expected = "Power-gated (\u2265 3 expected edited KO reads)", depth = "Depth-gated (\u2265 10 KO reads)")
V8_GATE_LAB2  <- c(expected = "Power-gated\n\u2265 3 expected edited KO reads", depth = "Depth-gated\n\u2265 10 KO reads")
V8_GATE_SHORT <- c(expected = "power-gated", depth = "depth-gated")
V8_TRUTH_LABELS <- setNames(c(
  sprintf("Ground truth, \u2265 %d%% editing reduction in KO, power-gated", c(25, 50, 75)),
  sprintf("Ground truth, \u2265 %d%% editing reduction in KO, depth-gated", c(25, 50, 75)),
  "REDIportal, \u2265 10 reads, \u2265 1 variant read", "REDIportal, \u2265 10 reads, \u2265 3 variant reads",
  "REDIportal, \u2265 10 reads, \u2265 3 variant reads, \u2265 5% editing level"), V8_TRUTH_IDS)
V8_TRUTH_COLOURS <- setNames(c("#99C7E0", "#0072B2", "#003F62", "#99D8C7", "#009E73", "#00573F",
                               "#F5D999", "#E69F00", "#7F5700"), V8_TRUTH_IDS)
V8_SUBTITLE   <- paste0("Call sets: sites near other mismatch types removed.\n",
                        "KO ground truths: sites near other mismatch types and PGP1 SNVs removed (no region mask). ",
                        "REDIportal: editing level \u2265 1%, PGP1 SNVs removed.\n",
                        "All truths: dbSNP common and multiallelic sites removed. ",
                        if (V8_SCORING == "pwcov")
                          paste0("Recall per unit: positives with \u2265 10 reads in that unit ",
                                 "(power-gated: also coverage \u00D7 NTC editing level \u2265 3).\n")
                        else "Recall per unit: positives with \u2265 10 reads in that unit.\n",
                        sprintf("Cell types: the %d with \u2265 %d truth sites under every truth", V8_N_CELLTYPES, V8_CT_MIN_SITES))
spearman_lab  <- function(rho, digits = 2) sprintf(paste0("Spearman's \u03C1 = %.", digits, "f"), rho)
