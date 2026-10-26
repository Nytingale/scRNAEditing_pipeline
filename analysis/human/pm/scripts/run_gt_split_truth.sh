#!/bin/bash
# =============================================================================
# run_gt_split_truth.sh  [submit|links|scan|cands|pile|eval|truth]   (2026-10-13)
#
# Split-by-cells test: the v8 truth (v7_nomm_pgp1) rebuilt from NTC cell half A only, with every
# rule frozen (FROZEN_RULES.tsv; nearby-mismatch w1 = 200, the full-data choice for both gates).
# Same chain as the full truth, read in place:
#   scan/cands/pile  as run_gt_spectrum12_v6.sh on bams/ = half-A NTC UMI BAMs (merged + T1-T4)
#                    and the full KO arm UMI BAMs (KO cells are not split)
#   eval             gt_spectrum12_eval_v6.R, both gates (parity warnings vs the full truth expected)
#   truth            gt_truth_agtc_test_v6.R (GT_TRUTH_ONLY) -> gt_truth_v7.R -> gt_truth_v7_variants.R
# Inputs : gt_eval_v8/split_cells/bams/umi/*_A.bam (split_gt_bams_halves.sh), 02_gt_bams_umi_v6/merged_ADAR*_KO.bam
# Outputs: gt_eval_v8/split_cells/{spectrum12,agtc_test,truth_v7,truth}/ ; truth/v7_nomm_pgp1/gt_v7_truth_sites.tsv.gz
#$ -cwd
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
SELF="$B/human/pm/scripts/run_gt_split_truth.sh"
PM="$B/human/pm"
S="$PM/gt_eval_v8/split_cells"
OUT="$S/spectrum12"
BAMDIR="$S/truth_bams"
FASTA="$B/human/ref_annots/hg38.analysisSet.fa"
REGIONS="$PM/gt_eval_v6/scan/regions.tsv"
SNPBED="$PM/gt_eval_v6/snp/dbSNP156_common.bed"
PARSER="$B/human/scripts/mismatch_site_count.pl"
UNITS=(merged_NTC_ctrl merged_ADAR1_KO merged_ADARB1_KO T1_NTC_ctrl T2_NTC_ctrl T3_NTC_ctrl T4_NTC_ctrl)
export GT_SP12_DIR="$OUT" GT_AGTC_DIR="$S/agtc_test" GT_V7_TRUTH_DIR="$S/truth_v7" GT_V7_VARIANTS_DIR="$S/truth"
export GT_TRUTH_ONLY=1 GT_W1_FIXED=200

mkdir -p "$OUT"/{scan,pileups,eval,eval_depth} "$BAMDIR"
cd "$B"
STAGE="${1:-submit}"
NREG=$(wc -l < "$REGIONS")

if [ "$STAGE" = "submit" ]; then
  HOLD=${HOLD:-}
  npile=$(( NREG * ${#UNITS[@]} ))
  jl=$(qsub -terse -N spl_links -q slice ${HOLD:+-hold_jid $HOLD} "$SELF" links)
  js=$(qsub -terse -N spl_scan -q slice -t 1-"$NREG" -hold_jid "$jl" "$SELF" scan | cut -d. -f1)
  jc=$(qsub -terse -N spl_cands -q slice -l h_vmem=32G -hold_jid "$js" "$SELF" cands)
  jp=$(qsub -terse -N spl_pile -q slice -t 1-"$npile" -hold_jid "$jc" "$SELF" pile | cut -d. -f1)
  je=$(qsub -terse -N spl_eval -q slice -l h_vmem=96G -hold_jid "$jp" "$SELF" eval)
  jt=$(qsub -terse -N spl_truth -q slice -l h_vmem=96G -hold_jid "$je" "$SELF" truth)
  echo "links $jl | scan $js (1-$NREG) | cands $jc | pile $jp (1-$npile) | eval $je | truth $jt"
  exit 0
fi

if [ "$STAGE" = "links" ]; then
  for u in merged_NTC_ctrl T1_NTC_ctrl T2_NTC_ctrl T3_NTC_ctrl T4_NTC_ctrl; do
    f="$S/bams/umi/${u}_A.bam"; [ -s "$f" ] && [ -s "$f.bai" ] || { echo "ERROR: missing $f(.bai)" >&2; exit 1; }
    ln -sf "$f" "$BAMDIR/$u.bam"; ln -sf "$f.bai" "$BAMDIR/$u.bam.bai"
  done
  for u in merged_ADAR1_KO merged_ADARB1_KO; do
    ln -sf "$PM/02_gt_bams_umi_v6/$u.bam" "$BAMDIR/$u.bam"; ln -sf "$PM/02_gt_bams_umi_v6/$u.bam.bai" "$BAMDIR/$u.bam.bai"
  done
  ls -lL "$BAMDIR"
  exit 0
fi

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
    | sort -u -k1,1 -k2,2n -S 8G > "$OUT/_positions_all.bed"
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
  T=$(mktemp -d "${TMPDIR:-/tmp}/spl.XXXXXX"); trap 'rm -rf "$T"' EXIT
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

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
