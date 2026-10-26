# plotting_functions.R
# Manuscript colour maps and shared ggplot2 helpers for scRNAEditing figures.
# Source this file from every script in plotting/. Do not put figure code here.
#
# PROJECT RULES (see .cursor/rules/plotting-colours.mdc):
#   1. Never put a red and a green on the same plot (red-green CVD, ~8% of men).
#      Build scales from `okabe_ito`; use blue-vs-orange for two-class
#      contrasts; use `scale_fill_manual_cvd()` / `scale_colour_manual_cvd()`,
#      which refuse an unsafe scale. `assert_cvd_safe()` checks any vector.
#   2. THIS FILE is the single home for every colour map. Figure scripts must
#      reference a `cmap_*` / `okabe_ito` / `accent_colour` object by name and
#      must not define colours inline.
#   3. Shorten large axis labels with `shorten_axis` (160000 -> "160K").
#   4. x and y axis text must be the same size, and the two axis titles the
#      same size as each other. `theme_manuscript()` enforces both.
#   5. `theme_manuscript()` = theme_classic for a single panel;
#      `theme_manuscript(facets = TRUE)` = theme_bw for facetted panels.
#
# Slim rebuild (2026-09-05) from misc_scripts/plotting_functions.R:
# colour maps, axis labels, theme, mismatch-table helpers only.
# Old HPC plotFilteringProcess* functions were not copied.
#
# Figures: SVG only (save_panel). Dense scatters use scattermore so the
# point layer is rasterized inside the SVG at dpi >= raster_dpi_min (900).
# Grids: grid_y = vertical bars, grid_x = horizontal bars,
# grid_both = scatters. Each panel writes <name>.sources.tsv and
# appends figure_source_map.tsv.

# Text is drawn with the Liberation Sans files in plotting/fonts/ only, so SVGs and PNGs are identical on any
# machine. Must be set before the first graphics device opens.
local({
  here <- NULL
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) { here <- fr$ofile; break }
  if (!is.null(here)) {
    conf <- file.path(dirname(normalizePath(here)), "fonts", "fonts.conf")
    if (file.exists(conf)) Sys.setenv(FONTCONFIG_FILE = conf)
  }
})

# ggplot2 (4.0) orders legends by a hash of their title, and for a legend without a title it looks `title` up
# outside the plot, finding graphics::title, whose serialization changes with the packages loaded. Untitled
# legends then swapped places between runs. With this object in scope the hash, and the order, are fixed.
title <- NULL

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(data.table)
})

# ── Palettes ────────────────────────────────────────────────────────────────
colorblind <- as.character(palette.colors(palette = "Okabe-Ito"))

bear <- c(
  "#faa818", "#41a30d", "#fbdf72", "#367d7d", "#d33502", "#6ebcbc", "#37526d",
  "#916848", "#f5b390", "#342739", "#bed678", "#a6d9ee", "#0d74b6",
  "#60824f", "#725ca5", "#e0598b"
)

palette_named <- function(keys, pal = bear) {
  keys <- unique(as.character(keys))
  out <- setNames(rep(pal, length.out = length(keys)), keys)
  warn_cvd_unsafe(out, "palette_named()")
  out
}

# ── Colour-vision-deficiency (CVD) safety ───────────────────────────────────
# PROJECT RULE: never put a red and a green on the same plot. Red-green CVD
# affects ~8% of men, who cannot separate the two. Build categorical scales
# from `okabe_ito` (CVD-safe by construction) and prefer blue-vs-orange for any
# two-class contrast. All colour maps belong in THIS file, never inline in a
# figure script.
okabe_ito <- c(
  black      = "#000000",
  orange     = "#E69F00",
  skyblue    = "#56B4E9",
  green      = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  purple     = "#CC79A7"
)

# Hue classification used by the checker below. Greys/near-neutrals and very
# dark or very light colours are exempt: they read as neutral, not as "red" or
# "green", so grey-vs-anything is always safe.
cvd_hue_class <- function(colours) {
  rgb <- grDevices::col2rgb(colours)
  hsv <- grDevices::rgb2hsv(rgb)
  hue <- hsv["h", ] * 360
  sat <- hsv["s", ]
  val <- hsv["v", ]
  out <- rep("neutral", length(colours))
  chromatic <- sat >= 0.25 & val >= 0.2 & val <= 0.98
  out[chromatic & (hue <= 25 | hue >= 335)] <- "red"
  out[chromatic & hue >= 80 & hue <= 175]   <- "green"
  out
}

# TRUE when a set of colours contains both a red and a green.
has_red_and_green <- function(colours) {
  cls <- cvd_hue_class(colours)
  any(cls == "red") && any(cls == "green")
}

warn_cvd_unsafe <- function(values, context = "this scale") {
  if (length(values) && has_red_and_green(values)) {
    cls <- cvd_hue_class(values)
    warning(
      "CVD: ", context, " puts red and green on the same plot ",
      "(red: ", paste(unique(values[cls == "red"]), collapse = ", "),
      " | green: ", paste(unique(values[cls == "green"]), collapse = ", "),
      "). Project rule: recolour one side, e.g. from okabe_ito.",
      call. = FALSE
    )
  }
  invisible(values)
}

# Hard stop, for figure scripts that should never ship an unsafe scale.
assert_cvd_safe <- function(values, context = "this scale") {
  if (length(values) && has_red_and_green(values)) {
    cls <- cvd_hue_class(values)
    stop("CVD: ", context, " puts red and green on the same plot (red: ",
         paste(unique(values[cls == "red"]), collapse = ", "), " | green: ",
         paste(unique(values[cls == "green"]), collapse = ", "),
         "). Project rule: recolour one side, e.g. from okabe_ito.",
         call. = FALSE)
  }
  invisible(values)
}

# Drop-in replacements for scale_*_manual() that enforce the rule.
scale_fill_manual_cvd <- function(values, ...) {
  assert_cvd_safe(values, "scale_fill_manual")
  scale_fill_manual(values = values, ...)
}

scale_colour_manual_cvd <- function(values, ...) {
  assert_cvd_safe(values, "scale_colour_manual")
  scale_colour_manual(values = values, ...)
}
scale_color_manual_cvd <- scale_colour_manual_cvd

# ── Two-class contrasts (blue vs orange: the safest CVD pair) ───────────────
# Reference vs edited allele, e.g. WASP removed-read composition.
cmap_allele <- c(
  "Reference (A)" = "#0072B2",
  "Edited (G)"    = "#E69F00",
  "ref"           = "#0072B2",
  "alt"           = "#E69F00"
)

# Read sets in the WASP reference-bias diagnostic.
cmap_readset <- c(
  "all reads"   = "grey60",
  "WASP-passed" = "#0072B2"
)

# Site/coverage attrition: what was measurable before WASP vs what survives it.
# Grey is the baseline being lost, so the ink follows the argument.
#
# The surviving bar is PURPLE, not the blue of `cmap_readset`/`cmap_allele`.
# The attrition panel is published stacked directly beneath the by-allele panel,
# and reusing blue there would read as "Evaluable sites == reference-allele
# reads" — two different quantities in two different units (sites vs reads).
# Different question, different colour.
cmap_attrition <- c(
  "Covered"   = "grey60",
  "Evaluable" = "#CC79A7"
)

# Whether a read overlaps a catalogued (REDIportal) variant, i.e. whether STAR
# WASP evaluates it at all.
cmap_variant_overlap <- c(
  "covers REDIportal site"       = "#0072B2",
  "covers no REDIportal site"    = "grey65"
)

# Whether a PGP1 germline SNV position is catalogued in dbSNP156 common.
cmap_dbsnp_overlap <- c(
  "In dbSNP156 common"     = "#0072B2",
  "Not in dbSNP156 common" = "grey65"
)

# Single accent colour for annotation lines, medians and reference lines.
# Vermillion, not pure red: safe because it is never paired with a green.
accent_colour <- "#D55E00"
neutral_colour <- "grey40"

# Filter outcome of a site in schematics (kept vs removed)
cmap_filter_keep <- c(Kept = "#0072B2", Removed = "#D55E00")

# KO ground-truth coverage gates
cmap_gate <- c("Power-gated" = "#0072B2", "Depth-gated" = "#E69F00")

# Full-length vs 3' platform comparison (Ding et al. PBMC; two-class contrast)
cmap_3prime_platform <- c("Smart-seq2" = "#0072B2", "10x 3' v3" = "#D55E00")
cmap_3prime_sites <- c("Sites with \u2265 10 reads" = "#009E73", "Edited sites" = "#CC79A7")

# Ground truth site classes
cmap_truth_class <- c("Positive" = "#0072B2", "Informative negative" = "#D55E00")

# ── Mismatch types ──────────────────────────────────────────────────────────
cmap_varType <- c(
  AC = "grey90", AG = "#D55E00", AT = "grey90",
  CA = "grey90", CG = "grey90", CT = "grey90",
  GA = "grey90",
  GC = "grey90", GT = "grey90", TA = "grey90",
  TC = "grey90", TG = "grey90"
)

# ── Datasets (display names and short codes) ────────────────────────────────
cmap_datasetname <- c(
  "Fly Brain"                       = "#D51F26",
  "Human Brain (10x V2)"            = "#2E006E",
  "Human Brain (10x V3)"            = "#53377A",
  "Human Lung"                      = "#987CD9",
  "Human Skin and Kidney"           = "#B99CFF",
  "Human Skin and Kidney (long)"    = "#DCB5FF",
  "Human PBMC (10x V2)"             = "#20244A",
  "Human PBMC (10x V3)"             = "#1926A6",
  "Human PBMC (Smart-seq2)"         = "#3D4FF2",
  "Human PBMC (InDrops)"            = "#6978FF",
  "Human PBMC (CEL-Seq2)"           = "#69A0FF",
  "Human PBMC (Drop-Seq)"           = "#69C3FF",
  "Human PBMC (Seq-Well)"           = "#57EBFF",
  "Mouse Cortex (10x V2)"           = "#FFB300",
  "Mouse Cortex (Smart-seq2)"       = "#FFCC00",
  "Mouse Cortex (DroNc-seq)"        = "#FFE100",
  "Mouse Cortex (sci-RNA-seq)"      = "#FFF07A",
  "Octopus Visual System"           = "chocolate4",
  "Octopus Brain"                   = "#B08A74",
  "Xenopus Embryo"                  = "#007D34",
  "Xenopus Animal Cap"              = "#7DD06F",
  "Zebrafish Embryo"                = "#F6768E"
)

cmap_dataset <- c(
  SA      = "#D51F26",
  ICV2    = "#2E006E",
  ICV3    = "#53377A",
  DS      = "#987CD9",
  RG      = "#B99CFF",
  JZL10x  = "#20244A",
  JZLXV3  = "#1926A6",
  JZLSS2  = "#3D4FF2",
  JZLIND  = "#6978FF",
  JZLCEL  = "#69A0FF",
  JZLDRO  = "#69C3FF",
  JZLSEQ  = "#57EBFF",
  mJZL10X = "#FFB300",
  mJZLSS2 = "#FFCC00",
  mJZLDNC = "#FFE100",
  mJZLSCI = "#FFF07A",
  CMN     = "chocolate4",
  ES      = "#B08A74",
  AMK     = "#007D34",
  KNN     = "#7DD06F",
  ACM     = "#F6768E",
  PM      = "#916848"
)

# ── Callers (pipeline slugs and manuscript labels) ──────────────────────────
cmap_method <- c(
  bcftools                               = "#00AAFF",
  reditools3                             = "#E3BE00",
  gatk                                   = "#C70021",
  jacusa2                                = "#6A3FA0",
  scallele                               = "#8C510A",
  cellsnp                                = "#F07800",
  "REDItools3"                           = "#E3BE00",
  "REDItools3 - Original Stranding"      = "#E3BE00",
  "REDItools3 - Stranded by annotations" = "#E3BE00",
  "REDItools3\n(Original stranding)"     = "#E3BE00",
  "REDItools3\n(Stranded by\nannotations)" = "#E3BE00",
  GATK                                   = "#C70021",
  JACUSA2                                = "#6A3FA0",
  scAllele                               = "#8C510A",
  "Cellsnp-lite"                         = "#F07800",
  BCFtools                               = "#00AAFF"
)

# ---- Editing-decrease bins (KO vs control) ---------------------------------
# Ordered levels + colour map for the "how far did editing drop in the KO" plots
# (human/pm/gt_eval/true_def_explore/plot_editing_decrease_*.R). Each evaluable
# A>G site is placed in one bin by the size of its editing decrease in the KO.
# Colour encodes severity: grey = NOT decreased (KO editing >= control), then a
# ColorBrewer Blues gradient that deepens as the decrease grows, up to the
# darkest blue = editing fully abolished in the KO (decreased to 0).
# `cmap_decrease_bin` keys MUST match `decrease_bin_levels` (the factor levels
# used in those scripts) exactly, embedded newlines included.
decrease_bin_levels <- c(
  "Not decreased\n(\u2265 ctrl)",   # KO editing >= control  -> grey
  "<5%", "5\u2013<10%", "10\u2013<15%",
  "15\u2013<20%", "20\u2013<25%", "\u226525%",
  "Decreased to 0"                  # KO editing == 0 at >=10x -> darkest blue
)
cmap_decrease_bin <- setNames(
  c("grey80", "#C6DBEF", "#9ECAE1", "#6BAED6",
    "#4292C6", "#2171B5", "#08519C", "#08306B"),
  decrease_bin_levels
)

# ---- Caller-name standardisation -------------------------------------------
# Canonical display names for the editing-site callers, so every figure/legend
# uses one spelling regardless of the raw slug (bcftools, BCFTOOLS, "jacusa-2",
# reditools_3, ...). Matching is case-insensitive and ignores separators.
# Usage:  df[, caller := std_caller(caller)]   or   labels = std_caller(levels)
caller_display <- c(
  bcftools   = "BCFtools",
  reditools3 = "REDItools3",
  reditools  = "REDItools",
  gatk       = "GATK",
  gatkhc     = "GATK",
  jacusa2    = "JACUSA2",
  jacusa     = "JACUSA2",
  scallele   = "scAllele",
  redml      = "RED-ML",
  sprint     = "SPRINT"
)

std_caller <- function(x) {
  raw <- as.character(x)
  key <- gsub("[^a-z0-9]", "", tolower(trimws(raw)))   # strip spaces/-/_/. etc
  out <- unname(caller_display[key])
  out[is.na(out)] <- raw[is.na(out)]                    # keep unknowns verbatim
  out
}

smap_method <- c(
  bcftools                               = 16,
  reditools3                             = 15,
  gatk                                   = 17,
  jacusa2                                = 18,
  scallele                               = 8,
  "REDItools3"                           = 15,
  "REDItools3 - Original stranding"      = 15,
  "REDItools3\n(Original stranding)"     = 15,
  "REDItools3 - Stranded by annotations" = 18,
  "REDItools3\n(Stranded by\nannotations)" = 18,
  GATK                                   = 17,
  JACUSA2                                = 18
)

# ── Resolution / BAM type ───────────────────────────────────────────────────
cmap_bulkType <- c(
  "Sample level"    = "#5E0800",
  "Cell type level" = "#D65347",
  "Per cell level"  = "#EDB5AF"
)

# Old manuscript labels (keep for reproducing submitted panels)
cmap_dedupType <- c(
  "Unprocessed"                         = "darkgreen",
  "De-duplicated"                       = "yellowgreen",
  "Recalibrated"                        = "darkorange1",
  "Recalibrated (Q20)"                  = "goldenrod1",
  "Annotated"                           = "steelblue4",
  "De-duplicated Annotated"             = "steelblue1",
  "Cell type pseudobulked"              = "steelblue4",
  "De-duplicated\ncell type pseudobulked" = "steelblue1"
)

# Dedup-audit panels (publication/teratoma_gt): blue / orange / skyblue —
# intentionally not the green cluster in cmap_bamType, so Input vs coordinate
# vs UMI stay separable under red-green CVD.
# Dedup-audit: grey for unprocessed baseline; warm hues for the two dedup
# strategies. NOT okabe-ito.
cmap_dedup_reads <- c(
  "Unprocessed"                  = "#AAAAAA",
  "De-duplicated (coordinate)"   = "#F28E2B",
  "De-duplicated (UMI)"          = "#4E79A7"
)
cmap_dedup_pct <- c(
  "MarkDuplicates" = "#F28E2B",
  "UMI"            = "#4E79A7"
)

# Pipeline BAM type names (cell filter × dedup)
cmap_bamType <- c(
  "sample" = "#AAAAAA",
  "Unprocessed" = "#AAAAAA",

  "sampleRmDup" = "#8C510A",
  "De-duplicated (coordinate)" = "#8C510A",
  "De-duplicated\n(coordinate)" = "#8C510A",

  "sampleUmiDup" = "#F28E2B",
  "De-duplicated (UMI)" = "#F28E2B",
  "De-duplicated\n(UMI)" = "#F28E2B",

  "goodCells"              = "#2E8B57",
  "goodCellsRmDup"         = "#66C2A5",
  "goodCellsUmiDup"        = "#8DD3C7",
  "annotated"              = "steelblue4",
  "annotatedRmDup"         = "steelblue1",
  "annotatedUmiDup"        = "#87CEEB",
  "celltype"               = "#37526D",
  "celltypeRmDup"          = "#6ebcbc",
  "celltypeUmiDup"         = "#a6d9ee",
  "perCell"                = "#C2185B",
  "perCellRmDup"           = "#E0598B",
  "perCellUmiDup"          = "#F4B4CE",

  "recal" = "#2171B5",
  "Recalibrated"  = "#2171B5",

  "recalq20" = "#6BAED6",
  "Recalibrated (Q20)" = "#6BAED6",
  "Recalibrated\n(Q20)" = "#6BAED6",

  "recaldefault" = "#08306B",
  "Recalibrated (default)" = "#08306B",
  "Recalibrated\n(default)" = "#08306B"
)

bam_type_levels <- names(cmap_bamType)

cmap_strandType <- c(
  "Stranded by\nannotations" = "darkorchid4",
  "Method's own\nstranding"  = "#EBCFFF"
)

cmap_strandMode <- c(
  annotation = "darkorchid4",
  hybrid     = "#EBCFFF"
)

cmap_qualType <- c("20" = "skyblue1", "30" = "navy")

cmap_geneRegion <- c(
  upstream              = "yellowgreen",
  downstream            = "limegreen",
  "upstream;downstream" = "seagreen",
  UTR5                  = "indianred1",
  "5UTR"                = "indianred1",
  UTR3                  = "orangered2",
  "3UTR"                = "orangered2",
  "UTR5;UTR3"           = "firebrick4",
  exonic                = "steelblue1",
  exons                 = "steelblue1",
  intronic              = "royalblue",
  introns               = "royalblue",
  splicing              = "dodgerblue4",
  ncRNA_exonic          = "gold",
  ncRNA_intronic        = "goldenrod1",
  ncRNA_splicing        = "darkgoldenrod3",
  intergenic            = "darkgrey"
)

cmap_subsamples <- c(
  "Subsample 1" = "#FAA818",
  "Subsample 2" = "#56B4E9",
  "Subsample 3" = "#009E73"
)

cmap_platform <- c(
  "10x V2"      = "#A6CDE2",
  "10x V3"      = "#688EC1",
  "10xv2"       = "#A6CDE2",
  "10xv3"       = "#688EC1",
  "InDrops"       = "#F59899",
  "Drop-Seq"    = "#FFD500",
  "CEL-Seq2"    = "#E11E26",
  "Smart-seq2"  = "#89288F",
  "smartseq2"     = "#89288F",
  "Seq-Well"    = "#B15928",
  "DroNc-seq"   = "#74C476",
  "sci-RNA-seq" = "#8C8C8C"
)

cmap_teratoma <- c(T1 = "#D55E00", T2 = "#0072B2", T3 = "#009E73", T4 = "#CC79A7",
                   ter1 = "#D55E00", ter2 = "#0072B2", ter3 = "#009E73", ter4 = "#CC79A7")

# Jitter-point palettes for the method-comparison panels (not okabe-ito).
# Teratoma: 4 hues from Tableau 10, all green-free.
cmap_teratoma_jitter <- c(
  T1 = "#4E79A7",  # steel blue
  T2 = "#F28E2B",  # orange
  T3 = "#B07AA1",  # mauve
  T4 = "#EDC948"   # gold
)

# Cell type: 23 hues grouped by lineage so biologically related types share a
# colour family. Green-free throughout (CVD-safe with the reds in the muscle
# group). Within each family, the broadest/parent type is the most saturated
# and derivatives are lighter or shifted.
#
#   Mesenchymal / fibroblast (warm orange-brown, 5 types)
#   Muscle                   (rose-red, 3 types)
#   Neural                   (blue-violet, 4 types)
#   Epithelial               (teal / slate, 4 types)
#   Neural crest             (pale blue / deep violet, 2 types)
#   Other                    (gold / grey / magenta / plum / lilac, 5 types)
cmap_teratoma_celltype <- c(
  # ── Mesenchymal / fibroblast ──
  MSC_Fib                 = "#B8651A",   # burnt orange (parent)
  Adipogenic_MSC_Fib      = "#E8943A",   # golden orange
  Cycling_MSC_Fib         = "#F2BA6B",   # light amber
  Chondrogenic_MSC_Fib    = "#C4813A",   # caramel
  MyoFib                  = "#7D4E24",   # dark brown
  # ── Muscle ──
  Cardiac_Skeletal_Muscle = "#8B1A2B",   # dark crimson
  Muscle_Prog             = "#D35070",   # medium rose
  Smooth_Muscle           = "#F0A0B0",   # light pink
  # ── Neural ──
  Radial_Glia             = "#2E4482",   # dark navy-blue
  Early_Neurons           = "#5B86D6",   # cornflower blue
  Retinal_Neurons         = "#7B52AB",   # medium purple
  Retinal_Epi             = "#B495D0",   # lavender
  # ── Epithelial ──
  Mid_Hindgut_Epi         = "#1B7A7A",   # dark teal
  Foregut_Epi             = "#6CC4B0",   # light seafoam
  Airway_Epi              = "#3FA3A3",   # mid teal
  Kidney_Prog             = "#2F4F5F",   # slate
  # ── Neural crest ──
  Schwann_Cells           = "#9DB9E8",   # pale blue
  Melanoblasts            = "#4A2D6E",   # deep violet
  # ── Other ──
  CycProg                 = "#E8C840",   # gold
  Pericytes               = "#9EA8AD",   # cool grey
  Immune                  = "#C44FC0",   # orchid / magenta
  HSC                     = "#7A1F72",   # plum
  Erythrocyte             = "#E9B3E4"    # pale lilac
)

# Guide-resolved teratoma CRISPR arms (Dailamy 2024 / GSE248941).
cmap_gt_arm <- c(
  ADAR1_KO   = "#D55E00",
  ADARB1_KO  = "#E69F00",
  ADARB2_KO  = "#56B4E9",
  AAVS1_ctrl = "#009E73",
  NTC_ctrl   = "#0072B2"
)

# ── Ground-truth definition (publication/gt_definition) ─────────────────────
# All four maps below are deliberately GREEN-FREE so they can sit beside the
# vermillion/orange used for the "not depleted" classes without tripping the
# red-green rule.

# v3 truth classes (GROUND_TRUTH_DEFINITION.md §2a.6). FP_unsplit is the label
# used when REDIportal membership is unavailable, so a missing catalogue is
# visible rather than silently collapsing every FP into "novel".
# Matched to cmap_gt_series: orange = depleted in KO, greys = not depleted.
# Class KEYS are the names fixed by GROUND_TRUTH_DEFINITION.md §2a.6; the
# figures relabel them for display (see GTDEF_gtdefinition_plot.R) because
# "FP_*" is a misnomer for a truth-set label.
cmap_gt_class <- c(
  "TRUE"                  = "#F1A340",
  "FP_known_not_depleted" = "#525252",
  "FP_novel"              = "#969696",
  "FP_unsplit"            = "#737373",
  "ambiguous"             = "grey85",
  "not_in_truth"          = "grey92"
)

# Aggregation tier. `union` is the combined site list the truth set is defined
# over, not a tier a caller can run at, hence the neutral grey.
# Pipeline token `dataset` is the arm-level merge (display label: "Arm").
cmap_gt_tier <- c(
  dataset  = "#0072B2",
  sample   = "#E69F00",
  celltype = "#CC79A7",
  union    = "grey40"
)

# Human-readable tier names for axis / facet / legend text.
label_gt_tier <- c(
  dataset  = "Arm",
  sample   = "Sample",
  celltype = "Cell type",
  union    = "Union"
)

# Figure / source_data basename suffix for a pipeline tier token.
# `dataset` → `arm` so filenames match the display label.
tier_file_stem <- function(tier) {
  t <- as.character(tier)
  t[t == "dataset"] <- "arm"
  t
}

cmap_tier_overlap <- c(
  shared         = "grey55",
  "arm only"     = "#0072B2",
  "sample only"  = "#E69F00"
)

# Step-4 alt-read rule, and the two consistency definitions compared in the
# ground-truth panels.
cmap_alt_rule <- c(count = "#0072B2", binom = "#E69F00")

cmap_consistency_mode <- c(strict = "#D55E00", evaluable = "#0072B2")

# Continuous scale for agreement/overlap heatmaps (Jaccard, containment).
# Grey-to-blue: sequential, green-free, and reads correctly in greyscale.
cmap_agreement <- c(low = "#F0F0F0", high = "#0072B2")

# Coverage strata, for panels where a claim depends on the depth cutoff chosen.
# Sequential light-to-dark blue so deeper strata read as "more resolved";
# green-free, and monotone in lightness so it survives greyscale printing.
cmap_cov_stratum <- c(
  ">= 1"   = "#9ECAE1",
  ">= 2"   = "#6BAED6",
  ">= 3"   = "#4292C6",
  ">= 5"   = "#2171B5",
  ">= 10"  = "#08519C",
  ">= 20"  = "#08306B",
  ">= 30"  = "#6BAED6",
  ">= 50"  = "#3182BD",
  ">= 100" = "#08519C",
  ">= 200" = "#08306B"
)

# Endpoints for sequential blue colourbars (e.g. the agreement heatmap), light
# through to very dark. Used on a FULL 0-1 scale: when the data all sit at the
# high end, every tile comes out dark, and that darkness is itself the message.
cmap_blue_ramp <- c(low = "#DEEBF7", high = "#08306B")

# ── Teratoma hues are RESERVED ──────────────────────────────────────────────
# PROJECT RULE: the four colours in `cmap_teratoma` (#D55E00, #0072B2, #009E73,
# #CC79A7) identify T1-T4 and must not be reused for any other category in the
# same figure. A reader who has just learnt "blue = T2" on one panel will read
# blue as T2 on the next one. The maps below therefore avoid all four, which is
# why they do not simply index into `okabe_ito`.

# Required depletion is ORDINAL, so it gets a sequential ramp rather than four
# categorical hues - which also sidesteps the teratoma palette entirely.
cmap_depletion <- c(
  "0.25" = "#C6DBEF",
  "0.5"  = "#6BAED6",
  "0.75" = "#2171B5",
  "0.9"  = "#08306B"
)

# Truth-set depletion classes.
#
# NAMING. Deliberately "depleted in KO", never "ADAR-dependent": a drop in
# editing in a knockout is an observed association, not proof the knockout
# caused it (off-target effects, clonal drift and coverage artefacts all
# produce the same signature). The panels should state what was measured.
# Equally never "false positive" at this level -- no method has made a call
# yet, so nothing can be wrong. FP only exists once a METHOD claims a site
# (gt_tier_scoring.tsv).
#
# COLOUR. The three depleted classes are NESTED subsets of increasing
# stringency, so they get a sequential warm ramp that darkens with stringency,
# and the non-depleted class gets a neutral dark grey. Green-free, so the reds
# are CVD-safe.
cmap_gt_series <- c(
  "Depleted in KO"                         = "#F1A340",
  "Depleted in KO, powered"                = "#D7301F",
  "Depleted in KO, CI excludes threshold"  = "#7F0000",
  "Not depleted in KO"                     = "#525252"
)

# ── GT v6 threshold sweeps (publication/gt_definition, v6 GTDEF scripts) ────
# PURPLE = overlap / agreement between site sets, everywhere in the GTDEF
# figures (sweep overlap lines and the consistency-agreement heatmap), so the
# hue keeps one meaning across a combined panel. Avoids the reserved teratoma
# hues and the orange/red/grey of cmap_gt_series.
cmap_gt_overlap = c(
  jaccard  = "#542788",   # Jaccard with the default TRUE set
  retained = "#B2ABD2"    # share of default TRUE sites retained
)
cmap_purple_ramp = c(low = "#F2F0F7", high = "#3F007D")

# Minimum-coverage floors, ordinal -> monotone light-to-dark blue (unlike
# cmap_cov_stratum, whose >= 30 / 50 entries restart the ramp).
cmap_cov_floor = c(
  ">= 5"   = "#9ECAE1",
  ">= 10"  = "#6BAED6",
  ">= 20"  = "#4292C6",
  ">= 50"  = "#2171B5",
  ">= 100" = "#08306B"
)

# Covered-A breakdown classes that exist only in the genome-wide v6 scan
# (extend cmap_gt_class; greys lighter than `ambiguous` = further from labelled).
cmap_gt_v6_extra = c(
  dbSNP_common_excluded = "#BDBDBD",
  below_candidate_floor = "#EFEFEF"
)

# GT v6 step figures (GTDEF_step*):
# TRUE-site COUNTS in 2-D sweep heatmaps: light-to-dark orange-red, the hue
# family of "Depleted in KO" in cmap_gt_series, so it reads as the same thing.
cmap_true_ramp = c(low = "#FEF0D9", high = "#B30000")
# Genome-wide POSITION counts (not truth sites) in heatmaps: neutral greys.
cmap_grey_ramp = c(low = "#F0F0F0", high = "#525252")
# % of TRUE sites in REDIportal, one line per sweep: dark teal, used for
# nothing else in the GTDEF figures.
cmap_gt_redi = c(rediportal = "#01665E")
# Minimum control edited reads 1..5 (ordinal) in the GT v6 joint diagnostic:
# light-to-dark green, a hue family not otherwise used in the GTDEF figures.
cmap_edit_reads_ramp = c("1" = "#C7E9C0", "2" = "#A1D99B", "3" = "#74C476", "4" = "#31A354", "5" = "#006D2C")
# Robustness (GTDEF_robust*): stability = overlap with neighbouring thresholds,
# so it reuses the purple overlap hue (dark = stable to 2-step moves).
cmap_gt_stable = c("Stable on all axes (1-step)" = "#3F007D", "Stable on continuous axes (1-step)" = "#9E9AC8",
                   "Not stable" = "#D9D9D9")
# GTDEF_step1_quality: pileup base/mapping-quality setting (dark = stricter).
cmap_bq_mq = c("BQ/MQ 30" = "#252525", "BQ/MQ 20" = "#A6A6A6")
# GTDEF_step2_universe: candidate sites in dbSNP (dbSNP156 common) vs not; used for nothing else.
cmap_dbsnp = c("Sites in dbSNP" = "#8C510A", "Sites not in dbSNP" = "#2166AC")
# Informative-KO explainer (GTDEF_robust_r5_informative_ko): merged-NTC editing
# level p of a site (ordinal, light = low editing); blues, used for nothing else.
cmap_ntc_edit_level = c("5%" = "#9ECAE1", "10%" = "#4292C6", "20%" = "#2171B5", "50%" = "#08306B")

# ── Strip labels placed on the Y side ───────────────────────────────────────
# PROJECT RULE: when facet strips sit on the y axis (facet_grid(..., switch =
# "y"), so the strip acts as a second axis title), drop the strip background and
# turn clipping off. The grey box reads as a panel of its own against the axis,
# and a two-line label gets clipped by the strip viewport without
# `strip.clip = "off"`.
#
#   p + facet_grid(metric ~ group, switch = "y") + theme_manuscript(facets = TRUE) +
#       strip_y_outside
strip_y_outside <- ggplot2::theme(
  strip.placement  = "outside",
  strip.background = ggplot2::element_blank(),
  strip.clip       = "off"
)

# Benchmark metrics. Three clearly SEPARATE hues: two shades of one colour make
# the reader group those two bars against the third, which is a grouping the
# data does not have.
cmap_gt_metric <- c(
  Precision = "#F1A340",
  Recall    = "#7B3294",
  F1        = "#08306B"
)

# ── Human-readable grid labels ──────────────────────────────────────────────
# PROJECT RULE: axis text must not be raw pipeline tokens. Nobody outside this
# project knows what "sampleRmDup" or "recalq20" means. Names follow the
# manuscript wording already used in `cmap_dedupType`, and the SNP wording
# follows nextflow/README.md's description of --snp_subtract_modes.
label_bamType <- c(
  sample           = "Unprocessed",
  sampleRmDup      = "De-duplicated (coordinate)",
  sampleUmiDup     = "De-duplicated (UMI)",
  recal            = "Recalibrated",
  recalq20         = "Recalibrated (Q20)",
  recaldefault     = "Recalibrated (default)",
  goodCells        = "QC-filtered cells",
  goodCellsRmDup   = "QC-filtered, de-duplicated (coordinate)",
  goodCellsUmiDup  = "QC-filtered, de-duplicated (UMI)",
  annotated        = "Annotated cells",
  annotatedRmDup   = "Annotated, de-duplicated (coordinate)",
  annotatedUmiDup  = "Annotated, de-duplicated (UMI)",
  celltype         = "Cell-type pseudobulk",
  celltypeRmDup    = "Cell-type pseudobulk, de-duplicated (coordinate)",
  celltypeUmiDup   = "Cell-type pseudobulk, de-duplicated (UMI)"
)

# PROJECT RULE: method names are written the way their authors write them,
# never the lowercase pipeline token.
label_caller <- c(
  bcftools   = "BCFtools",
  reditools3 = "REDItools3",
  gatk       = "GATK",
  jacusa2    = "JACUSA2",
  scallele   = "scAllele",
  cellsnp    = "Cellsnp-lite",
  redml      = "redML"
)
# Caller order for every plot and legend that is not ordered by rank.
CALLER_LEVELS <- c("BCFtools", "Cellsnp-lite", "GATK", "JACUSA2", "REDItools3", "scAllele")
caller_factor <- function(x, rev = FALSE) {
  y <- relabel(x, label_caller)
  lv <- c(CALLER_LEVELS, setdiff(sort(unique(y)), CALLER_LEVELS))
  factor(y, if (rev) rev(lv) else lv)
}

label_strandMode <- c(
  annotation = "Stranded by annotation",
  hybrid     = "Hybrid stranding",
  caller     = "Caller's own strand"
)

label_snpSubtract <- c(
  off = "dbSNP only",
  on  = "dbSNP + sample blacklist",
  na  = "dbSNP only"
)

# Apply a label map, passing through anything unmapped rather than turning it
# into NA -- a missing entry must be visible, not silently dropped.
relabel <- function(x, map) {
  x <- as.character(x)
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

# bam_type, for the ground-truth panels. Also avoids the reserved teratoma
# hues, so it is not built from `okabe_ito` by index. Green-free.
# BAM processing for ground-truth panels. "sample" (Input / Unprocessed) is
# deliberately grey: the unprocessed BAM is the baseline and grey signals
# "before anything was done". Other hues are from sequential blue (recal) and
# warm (dedup) ramps, NOT okabe-ito.
cmap_gt_bamtype <- c(
  recaldefault = "#08306B",
  recal        = "#2171B5",
  recalq20     = "#6BAED6",
  sampleUmiDup = "#F28E2B",
  sampleRmDup  = "#8C510A",
  sample       = "#AAAAAA"
)

# Each column gets its own categorical palette so the reader sees which
# preprocessing / caller / strand / SNP clusters at the top vs bottom.
cmap_proc <- c(
  "Unprocessed"                  = "#AFAFAF",
  "De-duplicated (coordinate)"   = "#487563",
  "De-duplicated\n(coordinate)"  = "#487563",
  "De-duplicated (UMI)"          = "#76A892",
  "De-duplicated\n(UMI)"         = "#76A892",
  "Recalibrated"                 = "#8055A6",
  "Recalibrated (Q20)"           = "#A47BBF",
  "Recalibrated (default)"       = "#C4A9D2"
)

cmap_strand <- c(
  "Stranded by annotation" = "#e89bc0",
  "Hybrid stranding"       = "#f5b673",
  "Caller's own strand"    = "#fff49e"
)

cmap_snp <- c(
  "dbSNP only"                    = "#CCCCCC",
  "dbSNP + sample blacklist"      = "#333333",
  "dbSNP + per-sample blacklist"  = "#333333",
  "dbSNP\n+ per-sample blacklist" = "#333333"
)

# Window filters (downstream steps 5-6): none = step 4; mismatch = step 5, sites with other
# mismatch types nearby (variable regions: artefact / SNP-prone) excluded; both = step 6,
# additionally only clustered sites retained
label_windowFilter <- c(
  none     = "None",
  mismatch = "Near other MM types excluded",
  both     = "Near other MM types excluded\n+ clustered"
)
cmap_window <- c(
  "None"                                        = "#DDEBF7",
  "Near other MM types excluded"                = "#6BAED6",
  "Near other MM types excluded\n+ clustered"   = "#08519C"
)
# legend text for figure captions (not drawn in figures)
caption_windowFilter <- paste0(
  "Local mismatch context: MM, mismatch. Near other MM types excluded = a site is removed when other mismatch types occur nearby;\n",
  "+ clustered = additionally, only sites with other candidate sites nearby are kept.")

# ── Axis / theme ────────────────────────────────────────────────────────────
# PROJECT RULE: shorten large axis labels. 160000 must read "160K", not
# "160,000". `shorten_axis` is the default for any count axis; axis_K / axis_M
# force a specific unit when a panel must stay comparable to a neighbour.
#   scale_y_continuous(labels = shorten_axis)
shorten_axis <- scales::label_number(scale_cut = scales::cut_short_scale())
axis_M <- function(l) paste0(round(l / 1e6, 1), "M")
axis_K <- function(l) paste0(round(l / 1e3, 1), "K")

# Same short-scale units, but keeping decimals. For IN-PLOT annotations, where
# `shorten_axis` rounds too hard and would print two different values
# identically (1,287,631 and 1,463,837 both become "1M"). Axis breaks are
# spaced out enough that `shorten_axis` is the right default for them.
shorten_num <- function(x, accuracy = 0.1) {
  scales::label_number(scale_cut = scales::cut_short_scale(),
                       accuracy = accuracy)(x)
}

my_theme <- theme(
  plot.title    = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  panel.grid    = element_blank()
)

# PROJECT RULE: theme_classic() for a single panel, theme_bw() for facets (the
# panel border and strip only earn their ink when there are panels to separate).
# PROJECT RULE: the x and y axes must use the SAME text size, and the same
# title size as each other. Before 2026-09-15 `my_theme` set axis.text.x = 12
# while y inherited base_size, so every figure had mismatched tick labels.
#
# PROJECT RULE: facet strip labels are ONE POINT LARGER than the axis text.
# A strip names the subset the whole panel is about, so it has to outrank the
# tick labels, but anything more turns the strip into a second title and
# competes with plot.title. Set here, not per figure, so every facetted panel
# in the manuscript agrees.
#
# PROJECT RULE: legend text matches the AXIS TEXT, and the legend title matches
# the AXIS TITLES. ggplot2's default is rel(0.8) of base_size, which makes the
# legend the smallest text on the panel — wrong, because a legend key is what
# decodes the colours and is read as often as a tick label. A legend entry and
# an axis tick label are the same kind of thing (both name a data value), so
# they get the same size.
#
#   theme_manuscript()              single panel  -> theme_classic
#   theme_manuscript(facets = TRUE) faceted       -> theme_bw
theme_manuscript <- function(base_size = 11, facets = FALSE,
                             axis_text_size = base_size,
                             axis_title_size = base_size + 1,
                             strip_text_size = axis_text_size + 1,
                             legend_text_size = axis_text_size,
                             legend_title_size = axis_title_size) {
  base <- if (isTRUE(facets)) {
    theme_bw(base_size = base_size)
  } else {
    theme_classic(base_size = base_size)
  }
  base +
    my_theme +
    theme(
      axis.text.x  = element_text(size = axis_text_size),
      axis.text.y  = element_text(size = axis_text_size),
      axis.title.x = element_text(size = axis_title_size),
      axis.title.y = element_text(size = axis_title_size),
      strip.text   = element_text(size = strip_text_size),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      legend.text  = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size),
      legend.position  = "bottom"
    )
}

# TRUE when a ggplot object is faceted, so a script can pass the right base
# theme without hardcoding it: p + theme_manuscript(facets = is_faceted(p))
is_faceted <- function(p) {
  f <- tryCatch(p$facet, error = function(e) NULL)
  !is.null(f) && !inherits(f, "FacetNull")
}

# Use after theme_manuscript(): vertical bars → grid_y, horizontal bars →
# grid_x, scatters → grid_both.
grid_y <- theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey85"),
  )

grid_x <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85"),
)

grid_both <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
)

# Default RNG for jitter and any other plotting randomness.
manuscript_seed <- 8787L
set.seed(manuscript_seed)

# Jitter first, then box_overlay() on top. Vertical categories: height = 0
# (jitter on x). Horizontal categories: width = 0 (jitter on y).
jitter_points <- function(orientation = c("vertical", "horizontal"),
                          seed = manuscript_seed, width = 0.2, height = 0.2,
                          ...) {
  orientation <- match.arg(orientation)
  if (orientation == "vertical") {
    geom_jitter(width = width, height = 0, seed = seed, ...)
  } else {
    geom_jitter(width = 0, height = height, seed = seed, ...)
  }
}

# Transparent fill, black outlines; hide outlier *glyphs* only
# (outlier.shape = NA). Do not set outliers = FALSE — that would drop
# those points from the box statistics and skew the hinges/whiskers.
# Unset mapped fill so a global aes(fill=...) cannot colour the boxes.
box_overlay <- function(width = 0.12, ...) {
  geom_boxplot(
    mapping = aes(fill = NULL),
    width = width, fill = NA, colour = "black", ...,
    outlier.shape = NA
  )
}

# ── Mismatch-count tables ───────────────────────────────────────────────────
all_vartypes <- function() {
  bases <- c("A", "C", "G", "T")
  CJ(ref = bases, alt = bases)[ref != alt, paste0(ref, alt)]
}

process_plotting_table <- function(df, setName) {
  dt <- as.data.table(df)
  dt[, label := paste0(
    format(n, big.mark = ","), "\n(",
    format(round(percent * 100, 2), nsmall = 2, trim = TRUE), "%)"
  )]
  dt[, setType := setName]
  dt[, colors := fifelse(varType == "AG", unname(cmap_varType[["AG"]]), "grey90")]
  if (nrow(dt) > 1) dt[varType != "AG", label := NA_character_]
  dt[]
}

add_missing_varTypes <- function(df) {
  dt <- as.data.table(df)
  missing <- setdiff(all_vartypes(), unique(as.character(dt$varType)))
  if (length(missing) == 0L) return(dt[order(varType)])
  extra <- data.table(varType = missing, n = 0L, percent = 0)
  rbind(dt, extra, fill = TRUE)[order(varType)]
}

collate_vartypes_tbl <- function(df, setTypeName) {
  dt <- as.data.table(df)
  if (!"varType" %in% names(dt) || nrow(dt) == 0L) {
    return(data.table(
      varType = all_vartypes(), n = 0L, percent = 0,
      totalSites = 0L, setType = setTypeName
    ))
  }
  tbl <- dt[, .(n = .N), by = varType]
  tbl <- add_missing_varTypes(tbl)
  tbl[is.na(n), n := 0L]
  tbl[, totalSites := sum(n)]
  tbl[, percent := ifelse(totalSites > 0, n / totalSites, 0)]
  tbl[, setType := setTypeName]
  tbl[]
}

# ─────────────────────────────────────────────────────────────────────────────
# >>> GRID COMPLETENESS — NO UNIQUE COMBINATION MAY BE SILENTLY DROPPED. <<<
#
# The benchmark is a grid over method x preprocessing x filter settings. If a
# panel does not map or facet EVERY varying dimension, the omitted one does not
# vanish politely — it corrupts the figure, and always silently. Measured on the
# 20-sample teratoma run (1280 rows, key unique), `snp_subtract` was absent from
# every panel, collapsing 1280 rows onto 740 plotted positions, in FOUR distinct
# failure modes:
#   geom_col + dodge  -> two bars at one position, one hides the other (p1)
#   geom_line(group)  -> run_label omitted it, so ONE line joined TWO different
#                        runs: a trajectory that does not exist (p2, p3)
#   merge()           -> non-unique key -> CARTESIAN product, 440 -> 880 points,
#                        half of them cross-paired nonsense (p4a, p4b)
#   geom_col(stack)   -> both runs STACKED, so bar heights are their SUM (p6)
#   geom_tile + text  -> two tiles per cell, one hidden, labels overlapping (p5)
#
# So every panel declares which dimensions it shows, and this warns if any
# varying dimension is unaccounted for. It WARNS rather than stops: it must never
# be the reason a multi-day pipeline run fails at the collation step.
#
# Candidates are the parameter columns of all_run_metrics.tsv. Constant ones are
# ignored, so adding a dimension to the grid automatically starts enforcing it
# instead of needing anyone to remember.
benchmark_grid_dims <- c("sample", "caller", "bam_type", "strand_mode",
                         "snp_subtract", "min_cov", "min_alt", "min_edit_level",
                         "jacusa2_min_score", "remove_100pct")

grid_dims_varying <- function(dt, candidates = benchmark_grid_dims) {
  candidates <- intersect(candidates, names(dt))
  candidates[vapply(candidates, function(cl) length(unique(dt[[cl]])) > 1L, logical(1))]
}

# >>> LEVEL LOSS IS THE OTHER WAY A COMBINATION DISAPPEARS. <<<
# `factor(x, levels = intersect(preferred, unique(x)))` is the idiom used to
# impose a house ordering, and it silently converts every value NOT in
# `preferred` to NA. Measured: `recaldefault` was absent from cmap_bamType, so
# 80 of 1280 runs became NA and were drawn as an unlabelled "NA" column in every
# panel. The dimension was present, so the dimension-level guard could not see
# it — hence this.
#
# Order by `preferred`, then APPEND anything else in sorted order. Unknown values
# are surfaced, never dropped, and a warning names them so the real fix (adding
# them to the colour map) still gets made.
factor_keep_all <- function(x, preferred, what = "value") {
  x <- as.character(x)
  present <- unique(x[!is.na(x)])
  extra <- sort(setdiff(present, preferred))
  if (length(extra)) {
    warning(sprintf(
      "%s(s) not in the house level order, appended rather than dropped: %s. Add them to the colour map.",
      what, paste(extra, collapse = ", ")), call. = FALSE)
  }
  factor(x, levels = c(intersect(preferred, present), extra))
}

# shown:      dimensions the panel distinguishes via aes/facet/group.
# aggregated: dimensions DELIBERATELY collapsed, e.g. sample when it is the
#   replicate dimension behind a median/IQR. Declaring one is the whole point —
#   an intended aggregation and an accidental overplot look identical in the
#   output, so the difference has to be stated in the code and is then auditable.
#   Anything neither shown nor declared is reported.
# position: optional: the columns that place a mark. If given, more than one row
#   per position means marks land on top of each other (or get stacked/summed).
assert_no_dropped_combos <- function(dt, shown, panel, position = NULL,
                                     aggregated = NULL,
                                     candidates = benchmark_grid_dims) {
  dims <- grid_dims_varying(dt, candidates)
  missing <- setdiff(dims, c(shown, aggregated))
  if (length(missing)) {
    warning(sprintf(
      paste0("panel '%s' does not distinguish varying grid dimension(s): %s. ",
             "Those combinations will be overplotted, fabricated into one group, ",
             "or summed. Map or facet them, or declare them in aggregated=."),
      panel, paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (length(aggregated)) {
    agg <- intersect(aggregated, names(dt))
    if (length(agg)) {
      message(sprintf("  %s: aggregating over %s (%s) — declared, not dropped",
                      panel, paste(agg, collapse = ", "),
                      paste(sprintf("n=%d", vapply(agg, function(cl)
                        length(unique(dt[[cl]])), integer(1))), collapse = ", ")))
    }
  }
  # NA in a dimension the panel plots means levels were lost, almost always by a
  # factor() whose levels= did not cover the data (see factor_keep_all).
  for (cl in intersect(shown, names(dt))) {
    n_na <- sum(is.na(dt[[cl]]))
    if (n_na) {
      warning(sprintf(
        paste0("panel '%s': %d of %d rows have NA '%s' — those combinations are ",
               "plotted unlabelled. A factor(levels=) almost certainly dropped a ",
               "level; use factor_keep_all()."),
        panel, n_na, nrow(dt), cl), call. = FALSE)
    }
  }
  if (!is.null(position)) {
    position <- intersect(position, names(dt))
    if (length(position)) {
      n_pos <- nrow(unique(dt[, position, with = FALSE]))
      if (n_pos < nrow(dt)) {
        warning(sprintf(
          "panel '%s': %d rows map onto %d distinct positions (%.1fx overplotted)",
          panel, nrow(dt), n_pos, nrow(dt) / n_pos), call. = FALSE)
      }
    }
  }
  invisible(setdiff(dims, shown))
}

# Rasterize point layers when a scatter has many points (keeps SVG small).
# Threshold is points on the plotted layer, not rows in a summarised table.
# Raster resolution is always dpi >= raster_dpi_min. Pass the same
# width/height (inches) that will be used in save_panel().
raster_point_min <- 500L
raster_dpi_min <- 900L

# >>> RASTER BUDGET — this is an OOM guard, not a quality setting. <<<
# scattermore allocates the whole pixels=c(w,h) raster (RGBA, plus intermediate
# copies), so raster memory grows with PANEL AREA x dpi^2 and is independent of
# the number of points. Panels whose width scales with the sample count therefore
# explode on wide datasets:
#     3 arms (Replogle)    w =  9 in -> 8,100 x 5,400   =  44 M px   fine
#    20 samples (teratoma) w = 60 in -> 54,000 x 5,400  = 292 M px   OOM
#    43 units (celltype)   w = 129 in -> 116,100 x 5,400 = 627 M px  worse
# That is what killed COLLATE_RESULTS in nf_pm_ps twice (jobs 1693635 at a 4 GB
# cap, 1705676 at 32 GB), both times on benchmark_strand_mode_ag_scatter.svg,
# and it is why the 3-arm Replogle run never hit it.
#
# dpi is reduced toward RASTER_DPI_FLOOR to keep the raster inside the budget.
# 300 dpi is the publication floor, and a 54,000-pixel-wide raster has no
# practical meaning anyway -- nothing prints or renders at that size. Panel
# geometry (width/height in inches) is UNCHANGED, so figure layout is identical;
# only the resolution of the embedded raster changes.
# >>> 600 dpi is a HARD FLOOR and outranks the budget. <<<
# Whenever a layer is rasterized it is written at >= raster_dpi_floor, 900
# preferred. The budget may only pull dpi DOWN TO the floor, never below it, so
# an oversized panel can never be silently published at a degraded resolution.
#
# 60 M px is chosen so that NO normal-sized panel is affected: the largest
# routine panel (8 x 8 in) is 7,200 x 7,200 = 51.8 M px and keeps the full
# 900 dpi, as does 9 x 6 (Replogle, 43.7 M px).
#
# CONSEQUENCE, deliberately left visible rather than papered over: a panel whose
# width scales with sample count can still exceed the budget at the floor
# (129 x 6 in at 600 dpi = 279 M px). The fix for THAT is the panel geometry, not
# the resolution -- `3 * n_samples` inches is wrong for facet_wrap, which wraps
# into a grid, so 20 facets ask for a 60-inch single-row strip. raster_over_budget()
# reports it so such a panel is caught rather than discovered by an OOM.
raster_dpi_floor <- 600L
raster_px_max    <- 60e6

# Effective dpi for a panel of this size: raster_dpi_min when it fits the budget,
# scaled down toward raster_dpi_floor when it does not.
raster_effective_dpi <- function(width_in, height_in, dpi = raster_dpi_min) {
  w <- as.numeric(width_in)[1];  h <- as.numeric(height_in)[1]
  if (is.na(w) || w <= 0) w <- 8
  if (is.na(h) || h <= 0) h <- 8
  dpi <- as.integer(dpi)[1]
  if (is.na(dpi)) dpi <- raster_dpi_min
  fit <- sqrt(raster_px_max / (w * h))
  as.integer(max(raster_dpi_floor, min(dpi, floor(fit))))
}

raster_pixels <- function(width_in, height_in, dpi = raster_dpi_min) {
  w <- as.numeric(width_in)[1]
  h <- as.numeric(height_in)[1]
  if (is.na(w) || w <= 0) w <- 8
  if (is.na(h) || h <= 0) h <- 8
  # NB no `max(dpi, raster_dpi_min)` clamp here any more: that is what made the
  # budget unenforceable, since it silently restored 900 dpi. The floor lives in
  # raster_effective_dpi() instead.
  dpi <- raster_effective_dpi(w, h, dpi)
  c(as.integer(ceiling(w * dpi)), as.integer(ceiling(h * dpi)))
}

# TRUE when a panel cannot fit the budget even at the dpi floor -- i.e. its
# GEOMETRY is the problem. Callers should shrink the panel, not the resolution.
raster_over_budget <- function(width_in, height_in) {
  w <- as.numeric(width_in)[1];  h <- as.numeric(height_in)[1]
  if (is.na(w) || w <= 0) w <- 8
  if (is.na(h) || h <= 0) h <- 8
  (w * raster_dpi_floor) * (h * raster_dpi_floor) > raster_px_max
}

points_layer <- function(n, mapping = NULL, pointsize = 2.2,
                         width = 8, height = 8, dpi = raster_dpi_min,
                         pixels = NULL, ...) {
  n <- as.integer(n)[1]
  if (is.na(n)) n <- 0L
  if (is.na(as.integer(dpi)[1])) dpi <- raster_dpi_min
  min_px <- raster_pixels(width, height, dpi)
  if (is.null(pixels)) {
    pixels <- min_px
  } else {
    pixels <- pmax(as.integer(pixels)[seq_len(2)], min_px)
  }
  # The budget binds even when pixels= was passed explicitly, otherwise a caller
  # could reintroduce the OOM -- but it may never push below the dpi floor.
  # Aspect ratio is preserved.
  floor_px <- c(as.integer(ceiling(as.numeric(width)[1]  * raster_dpi_floor)),
                as.integer(ceiling(as.numeric(height)[1] * raster_dpi_floor)))
  total <- prod(as.numeric(pixels))
  if (is.finite(total) && total > raster_px_max) {
    pixels <- as.integer(pmax(floor(pixels * sqrt(raster_px_max / total)), 1L))
    pixels <- pmax(pixels, floor_px)   # floor wins
  }
  if (n >= raster_point_min && raster_over_budget(width, height)) {
    warning(sprintf(
      paste0("rasterized panel is %.0f x %.0f in: %.0f M px even at the %d dpi floor ",
             "(budget %.0f M). Resolution is NOT being reduced below the floor -- ",
             "shrink the PANEL instead. `k * n_samples` inches is wrong for ",
             "facet_wrap, which wraps into a grid."),
      as.numeric(width)[1], as.numeric(height)[1],
      prod(as.numeric(floor_px)) / 1e6, raster_dpi_floor, raster_px_max / 1e6))
  }
  if (n >= raster_point_min && requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(
      mapping = mapping, pointsize = pointsize, pixels = pixels,
      interpolate = TRUE
    )
  } else {
    if (n >= raster_point_min) {
      warning("scattermore is not installed; writing a vector point layer. ",
              "Install scattermore and record the version in SOFTWARE_VERSIONS.txt.")
    }
    geom_point(mapping = mapping, ...)
  }
}

points_were_rasterized <- function(n) {
  as.integer(n)[1] >= raster_point_min &&
    requireNamespace("scattermore", quietly = TRUE)
}

# Paths under the project root (REPO_DIR, set by the *_common.R scripts) are recorded relative to it.
repo_relative <- function(p) {
  if (!exists("REPO_DIR", envir = globalenv())) return(p)
  root <- paste0(normalizePath(get("REPO_DIR", envir = globalenv()), mustWork = FALSE), "/")
  ifelse(startsWith(p, root), substring(p, nchar(root) + 1L), p)
}

# Force SVG. Record source TSVs and figure_source_map.tsv so each published
# panel can be traced to its source data.
#
# `sidecar_dir` is where those provenance/source-data TSVs go. Default NULL
# puts them beside the SVG, which is how the legacy teratoma_gt panels were
# written. Pass a directory to keep the figure folder SVG-only: a deliverable
# should not mix rendered panels with tables (see publication/README.md).
save_panel <- function(plot, file, width, height, sources,
                       data = NULL, grid = NA_character_,
                       rasterized = FALSE, script = NA_character_,
                       sidecar_dir = NULL) {
  if (missing(sources) || length(sources) == 0L) {
    stop("save_panel() requires sources= (paths of TSVs that produced this panel)")
  }
  sources <- unique(normalizePath(as.character(sources), mustWork = FALSE))
  sources_rel <- repo_relative(sources)
  file <- sub("\\.(pdf|png|svg)$", "", file, ignore.case = TRUE)
  svg_file <- paste0(file, ".svg")
  dir.create(dirname(svg_file), recursive = TRUE, showWarnings = FALSE)

  device <- if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite
  } else {
    grDevices::svg
  }
  ggsave(svg_file, plot, width = width, height = height, limitsize = FALSE,
         device = device, bg = "white")

  stem <- basename(file)
  side_dir <- if (is.null(sidecar_dir)) dirname(svg_file) else sidecar_dir
  dir.create(side_dir, recursive = TRUE, showWarnings = FALSE)

  src_dt <- data.table(
    plot_file   = basename(svg_file),
    source_file = sources_rel,
    source_exists = file.exists(sources)
  )
  fwrite(src_dt, file.path(side_dir, paste0(stem, ".sources.tsv")), sep = "\t")

  if (!is.null(data)) {
    fwrite(as.data.table(data),
           file.path(side_dir, paste0(stem, ".source_data.tsv")), sep = "\t")
  }

  map_file <- file.path(side_dir, "figure_source_map.tsv")
  row <- data.table(
    plot_file   = basename(svg_file),
    plot_path   = repo_relative(normalizePath(svg_file, mustWork = FALSE)),
    script      = script,
    grid        = as.character(grid)[1],
    rasterized  = isTRUE(rasterized),
    raster_dpi  = if (isTRUE(rasterized)) as.integer(raster_dpi_min) else NA_integer_,
    n_source_rows = if (is.null(data)) NA_integer_ else nrow(as.data.table(data)),
    source_files = paste(sources_rel, collapse = ";"),
    created_at  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  if (file.exists(map_file)) {
    old <- fread(map_file)
    old <- old[plot_file != row$plot_file]
    row <- rbind(old, row, fill = TRUE)
  }
  fwrite(row, map_file, sep = "\t")
  invisible(svg_file)
}

caller_script <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) return(normalizePath(f[1], mustWork = FALSE))
  NA_character_
}
