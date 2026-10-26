#!/bin/bash
# =============================================================================
# run_gt_editlevel_by_bamtype.sh
#   [submit|beds|pile|collect|correlate]
#
# BEFORE v5 ground truth: at each tier's discovery site list, quantify editing
# level on Unprocessed / coordinate-dedup / UMI-dedup for control + KO arms
# (NTC_ctrl, ADAR1_KO, ADARB1_KO). Correlate across bam types → decide which
# BAM processing to use for KO-depletion labels (v5 Phase 1).
#
# Site lists (not TRUE sites — those do not exist until after this decision):
#   true_def_v5_full/_site_lists/{sample,dataset,celltype}.bed
#   = unique sites from that tier's discovery table only (no cross-tier union).
#
# BAMs: results_gt_persample (only tree with all three processing states).
# Units: 3 tiers × 4 teratomas × 3 arms × 3 bam_types = 108 pileups.
#
# Do NOT run samtools mpileup on the login node.
# pe×vmem: keep h_vmem≤32G with pe≤2 so Illumina maps to slice.
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 2
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00

set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:/usr/local/bin:/usr/bin:/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true

B="${SCRNA_DIR:-.}"
ENG="$B/human/scripts"
PMS="$B/human/pm/scripts"
SELF="$PMS/run_gt_editlevel_by_bamtype.sh"
V5="$B/human/pm/gt_eval/true_def_v5_full"
SITELISTS="$V5/_site_lists"
PSROOT="$B/human/pm/results_gt_persample"
FASTA="$B/human/ref_annots/hg38.analysisSet.fa"
OUT="$B/human/pm/gt_eval/editlevel_by_bamtype_v5"
UNITS="$OUT/_units.txt"
FIG="$OUT/figures"

TIERS=(sample dataset celltype)
TERATOMAS=(T1 T2 T3 T4)
# Control + the two KO arms used in truth (§2a.5). AAVS1/ADARB2 excluded.
ARMS=(NTC_ctrl ADAR1_KO ADARB1_KO)
BAMTYPES=(unprocessed coord_dedup umi_dedup)

mkdir -p "$OUT/beds" "$OUT/pileups" "$OUT/frags" "$FIG"
cd "$B"
STAGE="${1:-submit}"

bam_path() {
  local tt="$1" arm="$2" bt="$3"
  local unit="${tt}_${arm}"
  case "$bt" in
    unprocessed) echo "${PSROOT}/${unit}/01_merged/${unit}_merged.bam" ;;
    coord_dedup) echo "${PSROOT}/${unit}/02_dedup/${unit}_rmDup.bam" ;;
    umi_dedup)   echo "${PSROOT}/${unit}/02_dedup/${unit}_umiDup.bam" ;;
  esac
}

# ------------------------------------------------------------------ submit ---
if [ "$STAGE" = "submit" ]; then
  # No hold on v5 classify — this chain *feeds* the truth decision.
  {
    for tier in "${TIERS[@]}"; do
      for tt in "${TERATOMAS[@]}"; do
        for arm in "${ARMS[@]}"; do
          for bt in "${BAMTYPES[@]}"; do
            echo -e "${tier}\t${tt}\t${arm}\t${bt}"
          done
        done
      done
    done
  } > "$UNITS"
  n=$(wc -l < "$UNITS")
  echo "units: $n (3 tiers × 4 teratomas × 3 arms × 3 bam_types)"
  echo "site lists: $SITELISTS/{sample,dataset,celltype}.bed"
  j1=$(qsub -terse -N el_beds -q slice -l h_vmem=8G "$SELF" beds | cut -d. -f1)
  echo "1 beds      -> $j1"
  j2=$(qsub -terse -N el_pile -q slice -hold_jid "$j1" -t 1-"$n" -l h_vmem=16G \
       "$SELF" pile | cut -d. -f1)
  echo "2 pile      -> $j2 (array 1-$n)"
  j3=$(qsub -terse -N el_coll -q slice -hold_jid "$j2" -l h_vmem=48G \
       "$SELF" collect | cut -d. -f1)
  echo "3 collect   -> $j3"
  j4=$(qsub -terse -N el_corr -q slice -hold_jid "$j3" -pe threaded 2 -l h_vmem=32G \
       "$SELF" correlate | cut -d. -f1)
  echo "4 correlate -> $j4"
  exit 0
fi

# -------------------------------------------------------------------- beds ---
# Symlink / copy discovery site-list BEDs into OUT/beds (one per tier).
if [ "$STAGE" = "beds" ]; then
  echo "=== $(date) === discovery site-list BEDs (pre-truth) ==="
  for t in "${TIERS[@]}"; do
    src="$SITELISTS/${t}.bed"
    [ -s "$src" ] || { echo "ERROR: missing $src" >&2; exit 1; }
    dest="$OUT/beds/${t}_sites.bed"
    cp -f "$src" "$dest"
    echo "  $t: $(wc -l < "$dest") sites -> $dest"
  done
  exit 0
fi

# -------------------------------------------------------------------- pile ---
if [ "$STAGE" = "pile" ]; then
  [ -s "$UNITS" ] || { echo "ERROR: missing $UNITS" >&2; exit 1; }
  TASK="${SGE_TASK_ID:?}"
  n=$(wc -l < "$UNITS")
  [ "$TASK" -ge 1 ] && [ "$TASK" -le "$n" ] || {
    echo "ERROR: SGE_TASK_ID=$TASK out of 1-$n" >&2; exit 1; }
  read -r TIER TT ARM BT < <(sed -n "${TASK}p" "$UNITS")
  BED="$OUT/beds/${TIER}_sites.bed"
  BAM=$(bam_path "$TT" "$ARM" "$BT")
  tag="${TIER}_${TT}_${ARM}_${BT}"
  OUTF="$OUT/pileups/pileup_${tag}.tsv"
  [ -s "$BED" ] || { echo "ERROR: missing $BED" >&2; exit 1; }
  [ -s "$BAM" ] || { echo "ERROR: missing $BAM" >&2; exit 1; }
  echo "=== task $TASK: tier=$TIER $TT $ARM $BT sites=$(wc -l < "$BED") ==="
  bash "$ENG/gt_site_pileup.sh" \
    --bam "$BAM" --bed "$BED" --fasta "$FASTA" --out "$OUTF" \
    --libtype forward --min-bq 30 --min-mq 30
  FRAG="$OUT/frags/frag_${tag}.tsv"
  TIER="$TIER" TT="$TT" ARM="$ARM" BT="$BT" OUTF="$OUTF" FRAG="$FRAG" Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))
p <- fread(Sys.getenv("OUTF"), sep = "\t", showProgress = FALSE)
p[, site := paste0(chr, ":", pos)]
out <- p[, .(tier = Sys.getenv("TIER"),
             teratoma = Sys.getenv("TT"),
             arm = Sys.getenv("ARM"),
             bam_type = Sys.getenv("BT"),
             site, cov = as.numeric(coverage),
             alt = as.numeric(variant_cov),
             edit = as.numeric(editing_level))]
fwrite(out, Sys.getenv("FRAG"), sep = "\t")
cat("frag rows:", nrow(out), "\n")
RS
  exit 0
fi

# ----------------------------------------------------------------- collect ---
if [ "$STAGE" = "collect" ]; then
  echo "=== $(date) === collect editing-level frags ==="
  n_exp=$(wc -l < "$UNITS")
  n_got=$(ls "$OUT"/frags/frag_*.tsv 2>/dev/null | wc -l)
  [ "$n_got" -eq "$n_exp" ] || {
    echo "ERROR: $n_got/$n_exp frags present" >&2
    ls "$OUT"/frags/frag_*.tsv 2>/dev/null | wc -l
    exit 1
  }
  OUT="$OUT" Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))
out <- Sys.getenv("OUT")
files <- list.files(file.path(out, "frags"), pattern = "^frag_.*\\.tsv$",
                    full.names = TRUE)
d <- rbindlist(lapply(files, fread), use.names = TRUE)
fwrite(d, file.path(out, "editlevel_long.tsv.gz"), sep = "\t")
cat("wrote editlevel_long.tsv.gz rows=", nrow(d),
    " sites=", uniqueN(d$site), "\n")
print(d[, .N, by = .(tier, arm, bam_type)][order(tier, arm, bam_type)])
RS
  exit 0
fi

# --------------------------------------------------------------- correlate ---
# Decision plots: are Unproc / coord / UMI editing levels interchangeable?
# Also: does KO/ctrl depletion ratio agree across bam types?
if [ "$STAGE" = "correlate" ]; then
  echo "=== $(date) === correlate editing levels + depletion across bam_types ==="
  OUT="$OUT" FIG="$FIG" Rscript - <<'RS'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
ggsave_svg <- function(file, plot, width, height, ...) {
  file <- sub("\\.(pdf|png)$", ".svg", file, ignore.case = TRUE)
  if (!grepl("\\.svg$", file, ignore.case = TRUE)) file <- paste0(file, ".svg")
  device <- if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite
  } else {
    grDevices::svg
  }
  ggplot2::ggsave(file, plot, width = width, height = height,
                  device = device, bg = "white", limitsize = FALSE, ...)
  invisible(file)
}

out <- Sys.getenv("OUT")
fig <- Sys.getenv("FIG")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
d <- fread(file.path(out, "editlevel_long.tsv.gz"), showProgress = FALSE)
min_cov <- 10L
pairs <- list(
  c("unprocessed", "coord_dedup"),
  c("unprocessed", "umi_dedup"),
  c("coord_dedup", "umi_dedup")
)
lab_bt <- c(unprocessed = "Unprocessed",
            coord_dedup = "Coord-dedup",
            umi_dedup   = "UMI-dedup")

# ---- per tier × arm × teratoma correlations --------------------------------
rows <- list()
for (tier_i in unique(d$tier)) {
  for (arm_i in unique(d$arm)) {
    for (tt_i in unique(d$teratoma)) {
      for (pr in pairs) {
        a <- d[tier == tier_i & arm == arm_i & teratoma == tt_i & bam_type == pr[1],
               .(site, edit_a = edit, cov_a = cov)]
        b <- d[tier == tier_i & arm == arm_i & teratoma == tt_i & bam_type == pr[2],
               .(site, edit_b = edit, cov_b = cov)]
        m <- merge(a, b, by = "site")
        m <- m[cov_a >= min_cov & cov_b >= min_cov]
        if (nrow(m) < 20L) next
        rows[[length(rows) + 1L]] <- data.table(
          tier = tier_i, arm = arm_i, teratoma = tt_i,
          bam_a = pr[1], bam_b = pr[2],
          n_sites = nrow(m),
          pearson  = cor(m$edit_a, m$edit_b, method = "pearson",
                         use = "complete.obs"),
          spearman = cor(m$edit_a, m$edit_b, method = "spearman",
                         use = "complete.obs"),
          mae = mean(abs(m$edit_a - m$edit_b), na.rm = TRUE),
          median_abs_delta = median(abs(m$edit_a - m$edit_b), na.rm = TRUE))
      }
    }
  }
}
corr <- rbindlist(rows)
fwrite(corr, file.path(out, "editlevel_correlation_by_teratoma.tsv"), sep = "\t")

# ---- pooled across teratomas -----------------------------------------------
pool_rows <- list()
for (tier_i in unique(d$tier)) {
  for (arm_i in unique(d$arm)) {
    for (pr in pairs) {
      a <- d[tier == tier_i & arm == arm_i & bam_type == pr[1] & cov >= min_cov,
             .(edit_a = mean(edit), n_tt = .N), by = site]
      b <- d[tier == tier_i & arm == arm_i & bam_type == pr[2] & cov >= min_cov,
             .(edit_b = mean(edit), n_tt_b = .N), by = site]
      m <- merge(a, b, by = "site")
      if (nrow(m) < 20L) next
      pool_rows[[length(pool_rows) + 1L]] <- data.table(
        tier = tier_i, arm = arm_i, bam_a = pr[1], bam_b = pr[2],
        n_sites = nrow(m),
        pearson  = cor(m$edit_a, m$edit_b, method = "pearson"),
        spearman = cor(m$edit_a, m$edit_b, method = "spearman"),
        mae = mean(abs(m$edit_a - m$edit_b)),
        median_abs_delta = median(abs(m$edit_a - m$edit_b)))
    }
  }
}
pool <- rbindlist(pool_rows)
fwrite(pool, file.path(out, "editlevel_correlation_pooled.tsv"), sep = "\t")

cat("\n=== median Spearman across teratomas (tier × arm × pair) ===\n")
print(corr[, .(median_spearman = round(median(spearman), 4),
               median_pearson  = round(median(pearson), 4),
               median_mae      = round(median(mae), 4),
               n = .N),
           by = .(tier, arm, pair = paste(bam_a, "vs", bam_b))][
             order(tier, arm, pair)])

cat("\n=== pooled ===\n")
print(pool[, .(tier, arm, pair = paste(bam_a, "vs", bam_b),
               n_sites, spearman = round(spearman, 4),
               pearson = round(pearson, 4), mae = round(mae, 4))][
                 order(tier, arm, pair)])

# ---- scatter: edit level Unproc vs dedup -----------------------------------
for (pr in pairs) {
  pts <- rbindlist(lapply(unique(d$tier), function(tier_i) {
    rbindlist(lapply(unique(d$arm), function(arm_i) {
      a <- d[tier == tier_i & arm == arm_i & bam_type == pr[1] & cov >= min_cov,
             .(edit_a = mean(edit)), by = site]
      b <- d[tier == tier_i & arm == arm_i & bam_type == pr[2] & cov >= min_cov,
             .(edit_b = mean(edit)), by = site]
      m <- merge(a, b, by = "site")
      if (!nrow(m)) return(NULL)
      m[, `:=`(tier = tier_i, arm = arm_i)]
      m
    }))
  }))
  rho <- pool[bam_a == pr[1] & bam_b == pr[2],
              .(tier, arm, lab = sprintf("rho==%.3f", spearman))]
  p <- ggplot(pts, aes(edit_a, edit_b)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_point(alpha = 0.12, size = 0.5, stroke = 0) +
    facet_grid(arm ~ tier) +
    geom_text(data = rho, aes(x = 0.05, y = 0.95, label = lab),
              parse = TRUE, hjust = 0, vjust = 1, size = 2.8,
              inherit.aes = FALSE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = paste0("editing level (", lab_bt[[pr[1]]], ")"),
         y = paste0("editing level (", lab_bt[[pr[2]]], ")"),
         title = paste0("Discovery-site editing: ", lab_bt[[pr[1]]],
                        " vs ", lab_bt[[pr[2]]]),
         subtitle = "Per-tier site lists; mean across teratomas; cov ≥ 10 both") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  fn <- file.path(fig, sprintf("editlevel_scatter_%s_vs_%s.svg", pr[1], pr[2]))
  cat("wrote", ggsave_svg(fn, p, width = 9, height = 7.5), "\n")
}

# ---- depletion ratio: KO_edit / NTC_edit per site × bam_type ---------------
# Decision-relevant: does bam choice change which sites look depleted?
dep_rows <- list()
for (tier_i in unique(d$tier)) {
  for (bt_i in unique(d$bam_type)) {
    for (ko_i in c("ADAR1_KO", "ADARB1_KO")) {
      ntc <- d[tier == tier_i & arm == "NTC_ctrl" & bam_type == bt_i & cov >= min_cov,
               .(edit_ntc = mean(edit), cov_ntc = mean(cov)), by = site]
      k <- d[tier == tier_i & arm == ko_i & bam_type == bt_i & cov >= min_cov,
             .(edit_ko = mean(edit), cov_ko = mean(cov)), by = site]
      m <- merge(ntc, k, by = "site")
      m <- m[edit_ntc > 0]
      if (!nrow(m)) next
      m[, `:=`(tier = tier_i, bam_type = bt_i, ko_arm = ko_i,
               ratio = edit_ko / edit_ntc,
               depleted_50 = edit_ko <= 0.5 * edit_ntc)]
      dep_rows[[length(dep_rows) + 1L]] <- m
    }
  }
}
dep <- rbindlist(dep_rows, use.names = TRUE, fill = TRUE)
fwrite(dep, file.path(out, "editlevel_depletion_ratios.tsv.gz"), sep = "\t")

# Concordance of depleted_50 call: Unproc vs each dedup
conc_rows <- list()
for (tier_i in unique(dep$tier)) {
  for (ko_i in unique(dep$ko_arm)) {
    for (pr in pairs) {
      a <- dep[tier == tier_i & ko_arm == ko_i & bam_type == pr[1],
               .(site, dep_a = depleted_50, ratio_a = ratio)]
      b <- dep[tier == tier_i & ko_arm == ko_i & bam_type == pr[2],
               .(site, dep_b = depleted_50, ratio_b = ratio)]
      m <- merge(a, b, by = "site")
      if (nrow(m) < 20L) next
      conc_rows[[length(conc_rows) + 1L]] <- data.table(
        tier = tier_i, ko_arm = ko_i, bam_a = pr[1], bam_b = pr[2],
        n_sites = nrow(m),
        agree = mean(m$dep_a == m$dep_b),
        both_dep = mean(m$dep_a & m$dep_b),
        only_a = mean(m$dep_a & !m$dep_b),
        only_b = mean(!m$dep_a & m$dep_b),
        spearman_ratio = cor(m$ratio_a, m$ratio_b, method = "spearman",
                             use = "complete.obs"))
    }
  }
}
conc <- rbindlist(conc_rows)
fwrite(conc, file.path(out, "editlevel_depletion_concordance.tsv"), sep = "\t")
cat("\n=== depletion-call concordance (50% threshold) ===\n")
print(conc[, .(tier, ko_arm,
               pair = paste(bam_a, "vs", bam_b),
               n_sites, agree = round(agree, 3),
               spearman_ratio = round(spearman_ratio, 3))][
                 order(tier, ko_arm, pair)])

# Scatter of depletion ratios Unproc vs coord / umi
for (pr in list(c("unprocessed", "coord_dedup"),
                c("unprocessed", "umi_dedup"))) {
  pts <- rbindlist(lapply(unique(dep$tier), function(tier_i) {
    rbindlist(lapply(unique(dep$ko_arm), function(ko_i) {
      a <- dep[tier == tier_i & ko_arm == ko_i & bam_type == pr[1],
               .(site, ratio_a = ratio)]
      b <- dep[tier == tier_i & ko_arm == ko_i & bam_type == pr[2],
               .(site, ratio_b = ratio)]
      m <- merge(a, b, by = "site")
      if (!nrow(m)) return(NULL)
      m[, `:=`(tier = tier_i, ko_arm = ko_i)]
      m
    }))
  }))
  rho <- conc[bam_a == pr[1] & bam_b == pr[2],
              .(tier, ko_arm = ko_arm,
                lab = sprintf("rho==%.3f\nagree==%.2f",
                              spearman_ratio, agree))]
  p <- ggplot(pts, aes(ratio_a, ratio_b)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.3) +
    geom_point(alpha = 0.12, size = 0.5, stroke = 0) +
    facet_grid(ko_arm ~ tier) +
    geom_text(data = rho, aes(x = 0.02, y = 1.8, label = lab),
              hjust = 0, vjust = 1, size = 2.5, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
    labs(x = paste0("KO/NTC edit ratio (", lab_bt[[pr[1]]], ")"),
         y = paste0("KO/NTC edit ratio (", lab_bt[[pr[2]]], ")"),
         title = paste0("Depletion ratio: ", lab_bt[[pr[1]]],
                        " vs ", lab_bt[[pr[2]]]),
         subtitle = "Dashed lines = 50% depletion threshold; discovery site lists") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  fn <- file.path(fig, sprintf("depletion_ratio_%s_vs_%s.svg", pr[1], pr[2]))
  cat("wrote", ggsave_svg(fn, p, width = 9, height = 6.5), "\n")
}

# Summary bar: median Spearman edit-level Unproc vs dedup (NTC only — ctrl axis)
summ <- pool[arm == "NTC_ctrl" & bam_a == "unprocessed",
             .(tier, bam_b, spearman, mae, n_sites)]
p_sum <- ggplot(summ, aes(tier, spearman, fill = bam_b)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c(coord_dedup = "#4E79A7", umi_dedup = "#F28E2B"),
                    labels = c(coord_dedup = "vs Coord-dedup",
                               umi_dedup = "vs UMI-dedup")) +
  labs(x = NULL, y = "Spearman (editing level)",
       fill = NULL,
       title = "NTC editing-level agreement: Unprocessed vs dedup",
       subtitle = "Pooled across teratomas; cov ≥ 10; dashed = 0.95") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
cat("wrote", ggsave_svg(file.path(fig, "editlevel_agreement_summary.svg"), p_sum, width = 6.5, height = 4), "\n")
cat("done\n")
RS
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
