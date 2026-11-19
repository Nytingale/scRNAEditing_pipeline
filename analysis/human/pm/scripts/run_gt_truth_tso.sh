#!/bin/bash
# =============================================================================
# run_gt_truth_tso.sh  [submit|setup|scan|cands|pile|eval|truth|copy|robust]   (2026-10-29)
#
# TSO-clipped teratoma: the v8 truth (v7_nomm_pgp1) rebuilt from the clipped UMI BAMs (02_gt_bams_umi_v6), then the
# robustness sweeps on it. Same chain as run_gt_split_truth.sh, at the default (full-data) locations, with the w1
# search on (no GT_W1_FIXED):
#   setup            gt_eval_v6/: FROZEN_RULES.tsv, scan/regions.tsv, snp/dbSNP156_common.bed and the parity
#                    references (old v6 truth, sweep_infneg{,_depth} classes) linked from human/pm_noclip
#   scan/cands/pile  run_gt_spectrum12_v6.sh stages (merged NTC scan; pileups of the 7 truth units)
#   eval             gt_spectrum12_eval_v6.R, both gates (parity lines compare with the pre-clipping truth)
#   truth            gt_truth_agtc_test_v6.R (GT_TRUTH_ONLY) -> gt_truth_v7.R -> gt_truth_v7_variants.R
#   copy             v7_nomm_pgp1 -> gt_eval_v8/truth/ (+ support/), as gt_eval_v8/README.md
#   robust           gt_robustness_v8.R -> gt_eval_v8/robustness/
# Usage: bash human/pm/scripts/run_gt_truth_tso.sh   (login node; setup runs here, the rest is queued)
# =============================================================================
#$ -cwd
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
set -eo pipefail
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"

B="${SCRNA_DIR:-.}"
SELF="$B/human/pm/scripts/run_gt_truth_tso.sh"
PM="$B/human/pm"
NC="$B/human/pm_noclip"
V6="$PM/gt_eval_v6"
OUT="${GT_SP12_DIR:-$V6/spectrum12}"   # seed replicates: GT_SP12_DIR / GT_UMI_BAMDIR (umi_seed_reps.sh)
BAMDIR="${GT_UMI_BAMDIR:-$PM/02_gt_bams_umi_v6}"
FASTA="$B/human/ref_annots/hg38.analysisSet.fa"
REGIONS="$V6/scan/regions.tsv"
SNPBED="$V6/snp/dbSNP156_common.bed"
PARSER="$B/human/scripts/mismatch_site_count.pl"
UNITS=(merged_NTC_ctrl merged_ADAR1_KO merged_ADARB1_KO T1_NTC_ctrl T2_NTC_ctrl T3_NTC_ctrl T4_NTC_ctrl)
export GT_TRUTH_ONLY=1 TMPDIR="$B/tmp"
cd "$B"
STAGE="${1:-submit}"

if [ "$STAGE" = "submit" ]; then
  bash "$SELF" setup
  NREG=$(wc -l < "$REGIONS"); npile=$(( NREG * ${#UNITS[@]} ))
  js=$(qsub -terse -N tsot_scan -t 1-"$NREG" "$SELF" scan | cut -d. -f1)
  jc=$(qsub -terse -N tsot_cands -l h_vmem=32G -hold_jid "$js" "$SELF" cands)
  jp=$(qsub -terse -N tsot_pile -t 1-"$npile" -hold_jid "$jc" "$SELF" pile | cut -d. -f1)
  je=$(qsub -terse -N tsot_eval -l h_vmem=96G -hold_jid "$jp" "$SELF" eval)
  jt=$(qsub -terse -N tsot_truth -l h_vmem=96G -hold_jid "$je" "$SELF" truth)
  jy=$(qsub -terse -N tsot_copy -l h_vmem=4G -hold_jid "$jt" "$SELF" copy)
  jr=$(qsub -terse -N tsot_robust -l h_vmem=32G -hold_jid "$jy" "$SELF" robust)
  echo "scan $js (1-$NREG) | cands $jc | pile $jp (1-$npile) | eval $je | truth $jt | copy $jy | robust $jr"
  exit 0
fi

if [ "$STAGE" = "setup" ]; then
  mkdir -p "$V6"/{scan,snp,truth,scoring/sweep_infneg,scoring/sweep_infneg_depth} "$OUT"/{scan,pileups,eval,eval_depth}
  # Inputs: FROZEN_RULES.tsv (in the repository), regions + dbSNP BED (run_gt_truth_rebuild_v6.sh init / prep);
  # parity references (optional): the pre-clipping v6 truth and sweep classes. Linked from human/pm_noclip when there.
  [ -e "$V6/FROZEN_RULES.tsv" ] || cp "$NC/gt_eval_v6/FROZEN_RULES.tsv" "$V6/"
  for f in scan/regions.tsv snp/dbSNP156_common.bed truth/gt_truth_sites.tsv.gz \
           scoring/sweep_infneg/truth_sweep_classes.tsv.gz scoring/sweep_infneg_depth/truth_sweep_classes.tsv.gz; do
    [ -e "$V6/$f" ] || [ ! -e "$NC/gt_eval_v6/$f" ] || ln -s "$NC/gt_eval_v6/$f" "$V6/$f"
  done
  [ -s "$REGIONS" ] && [ -s "$SNPBED" ] || { echo "ERROR: missing $REGIONS or $SNPBED; build them first:" \
    "bash human/pm/scripts/run_gt_truth_rebuild_v6.sh init; qsub human/pm/scripts/run_gt_truth_rebuild_v6.sh prep" >&2; exit 1; }
  for u in "${UNITS[@]}"; do [ -s "$BAMDIR/$u.bam" ] && [ -s "$BAMDIR/$u.bam.bai" ] || { echo "ERROR: missing $BAMDIR/$u.bam(.bai)" >&2; exit 1; }; done
  exit 0
fi

NREG=$(wc -l < "$REGIONS")

if [ "$STAGE" = "scan" ]; then
  read -r CHUNK CHR S0 E REG <<< "$(sed -n "${SGE_TASK_ID:?}p" "$REGIONS")"
  rm -f "$OUT/scan/_ok_${CHUNK}"
  samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -f "$FASTA" "$BAMDIR/merged_NTC_ctrl.bam" \
    | perl "$PARSER" scan 5 1 0.01 | gzip -c > "$OUT/scan/${CHUNK}.cand.tsv.gz"
  touch "$OUT/scan/_ok_${CHUNK}"
  exit 0
fi

if [ "$STAGE" = "cands" ]; then
  nok=$(ls "$OUT"/scan/_ok_* 2>/dev/null | wc -l)
  [ "$nok" -eq "$NREG" ] || { echo "ERROR: $nok/$NREG scan chunks OK" >&2; exit 1; }
  cut -f1 "$REGIONS" | while read -r c; do zcat "$OUT/scan/${c}.cand.tsv.gz" | tail -n +2; done \
    | gzip -c > "$OUT/candidates_all.tsv.gz"
  zcat "$OUT/candidates_all.tsv.gz" | awk -F'\t' -v OFS='\t' '{print $1, $2 - 1, $2}' \
    | sort -u -k1,1 -k2,2n -S 8G -T "$TMPDIR" > "$OUT/_positions_all.bed"
  bedtools intersect -v -a "$OUT/_positions_all.bed" -b "$SNPBED" > "$OUT/positions_noSNP.bed"
  {
    printf "metric\tn\n"
    printf "candidate_rows\t%s\n" "$(zcat "$OUT/candidates_all.tsv.gz" | wc -l)"
    printf "positions_all\t%s\n" "$(wc -l < "$OUT/_positions_all.bed")"
    printf "positions_noSNP\t%s\n" "$(wc -l < "$OUT/positions_noSNP.bed")"
  } | tee "$OUT/cands_stats.tsv"
  exit 0
fi

if [ "$STAGE" = "pile" ]; then
  i=$(( SGE_TASK_ID - 1 )); u=$(( i / NREG )); r=$(( i % NREG + 1 ))
  UNIT="${UNITS[$u]}"
  read -r CHUNK CHR S0 E REG <<< "$(sed -n "${r}p" "$REGIONS")"
  F="$OUT/pileups/${UNIT}__${CHUNK}.tsv.gz"
  rm -f "$OUT/pileups/_ok_${UNIT}__${CHUNK}"
  T=$(mktemp -d "$TMPDIR/tsot.XXXXXX"); trap 'rm -rf "$T"' EXIT
  awk -v c="$CHR" -v s="$S0" -v e="$E" '$1 == c && $3 >= s && $3 <= e' "$OUT/positions_noSNP.bed" > "$T/pos.bed"
  if [ -s "$T/pos.bed" ]; then
    samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -l "$T/pos.bed" -f "$FASTA" "$BAMDIR/${UNIT}.bam" \
      | perl "$PARSER" site | gzip -c > "$F"
  else
    perl "$PARSER" site < /dev/null | gzip -c > "$F"
  fi
  touch "$OUT/pileups/_ok_${UNIT}__${CHUNK}"
  exit 0
fi

if [ "$STAGE" = "eval" ]; then
  n=$(ls "$OUT"/pileups/_ok_* 2>/dev/null | wc -l); nexp=$(( NREG * ${#UNITS[@]} ))
  [ "$n" -eq "$nexp" ] || { echo "ERROR: $n/$nexp pileup chunks OK" >&2; exit 1; }
  GT_COVER_GATE=expected Rscript "$B/human/scripts/gt_spectrum12_eval_v6.R"
  GT_COVER_GATE=depth Rscript "$B/human/scripts/gt_spectrum12_eval_v6.R"
  exit 0
fi

if [ "$STAGE" = "truth" ]; then
  Rscript "$B/human/scripts/gt_truth_agtc_test_v6.R"
  Rscript "$B/human/scripts/gt_truth_v7.R"
  Rscript "$B/human/scripts/gt_truth_v7_variants.R"
  exit 0
fi

if [ "$STAGE" = "copy" ]; then
  V=$PM/gt_eval_v7/variants/truth; D=$PM/gt_eval_v8/truth
  mkdir -p "$D/support"
  cp "$V/v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz" "$D/"
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "variant") c = i; print; next } $c == "v7_nomm_pgp1"' \
    "$V/truth_variant_counts.tsv" > "$D/gt_v8_truth_counts.tsv"
  cp "$V"/{truth_variant_counts,spectrum12_cleaning,nomm_w1_grid}.tsv "$V"/truth{,12}_site_flags.tsv.gz "$D/support/"
  cat "$D/gt_v8_truth_counts.tsv"
  awk -F'\t' 'NR == 1 || $NF == "TRUE" || /TRUE/' "$D/support/nomm_w1_grid.tsv" | head -20
  exit 0
fi

if [ "$STAGE" = "robust" ]; then
  Rscript "$B/human/scripts/gt_robustness_v8.R"
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
