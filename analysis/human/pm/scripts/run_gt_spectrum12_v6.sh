#!/bin/bash
# =============================================================================
# run_gt_spectrum12_v6.sh  [submit|scan|cands|pile|eval|plot]
#
# Mismatch spectrum of the v6 ground truth: the v6 truth rules (frozen, see
# gt_eval_v6/FROZEN_RULES.tsv) applied to ALL 12 mismatch types instead of A>G
# only, genome-wide and caller-independent, for every depletion threshold d of
# the step-6 sweep (25-90% in 5% steps). Type AG at d = 50% must reproduce the
# 1,692 v6 ground-truth sites (checked in `eval`).
#
#   scan   merged NTC UMI BAM, per 40 Mb chunk (same regions as v6 scan):
#          candidates per (site, strand, alt) at the v6 floor (cov >= 5,
#          alt >= 1, >= 1%) via human/scripts/mismatch_site_count.pl scan
#   cands  merge; drop dbSNP common positions (same BED as v6); union BED
#   pile   raw per-strand counts at the union, BQ/MQ 30, for the 7 units the
#          rules use (merged NTC / ADAR1_KO / ADARB1_KO, T1-T4 NTC), per chunk
#   eval   human/scripts/gt_spectrum12_eval_v6.R
#   plot   publication/gt_definition/scripts/GTDEF_spectrum12_plot.R
#
# Inputs (read in place): human/pm/02_gt_bams_umi_v6/*.bam, gt_eval_v6/scan/regions.tsv,
#   gt_eval_v6/snp/dbSNP156_common.bed, gt_eval_v6/truth/gt_truth_sites.tsv.gz
# Outputs: human/pm/gt_eval_v6/spectrum12/   Logs: logs/
# Do NOT run samtools on the login node; `submit` only calls qsub.
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00

set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"

B="${SCRNA_DIR:-.}"
SELF="$B/human/pm/scripts/run_gt_spectrum12_v6.sh"
V6="$B/human/pm/gt_eval_v6"
OUT="$V6/spectrum12"
BAMDIR="$B/human/pm/02_gt_bams_umi_v6"
FASTA="$B/human/ref_annots/hg38.analysisSet.fa"
REGIONS="$V6/scan/regions.tsv"
SNPBED="$V6/snp/dbSNP156_common.bed"
PARSER="$B/human/scripts/mismatch_site_count.pl"
UNITS=(merged_NTC_ctrl merged_ADAR1_KO merged_ADARB1_KO T1_NTC_ctrl T2_NTC_ctrl T3_NTC_ctrl T4_NTC_ctrl)

mkdir -p "$OUT"/{scan,pileups,eval}
cd "$B"
STAGE="${1:-submit}"
NREG=$(wc -l < "$REGIONS")

if [ "$STAGE" = "submit" ]; then
  npile=$(( NREG * ${#UNITS[@]} ))
  js=$(qsub -terse -N sp12_scan -q slice -t 1-"$NREG" "$SELF" scan | cut -d. -f1)
  jc=$(qsub -terse -N sp12_cands -q slice -l h_vmem=32G -hold_jid "$js" "$SELF" cands)
  jp=$(qsub -terse -N sp12_pile -q slice -t 1-"$npile" -hold_jid "$jc" "$SELF" pile | cut -d. -f1)
  je=$(qsub -terse -N sp12_eval -q slice -l h_vmem=96G -hold_jid "$jp" "$SELF" eval)
  jf=$(qsub -terse -N sp12_plot -q slice -l h_vmem=16G -hold_jid "$je" "$SELF" plot)
  echo "scan $js (1-$NREG) | cands $jc | pile $jp (1-$npile) | eval $je | plot $jf"
  exit 0
fi

if [ "$STAGE" = "scan" ]; then
  read -r CHUNK CHR S E REG <<< "$(sed -n "${SGE_TASK_ID:?}p" "$REGIONS")"
  rm -f "$OUT/scan/_ok_${CHUNK}"
  samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -f "$FASTA" "$BAMDIR/merged_NTC_ctrl.bam" \
    | perl "$PARSER" scan 5 1 0.01 | gzip -c > "$OUT/scan/${CHUNK}.cand.tsv.gz"
  touch "$OUT/scan/_ok_${CHUNK}"
  echo "$CHUNK: $(zcat "$OUT/scan/${CHUNK}.cand.tsv.gz" | tail -n +2 | wc -l) candidate (site, strand, alt) rows"
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
  read -r CHUNK CHR S E REG <<< "$(sed -n "${r}p" "$REGIONS")"
  F="$OUT/pileups/${UNIT}__${CHUNK}.tsv.gz"
  rm -f "$OUT/pileups/_ok_${UNIT}__${CHUNK}"
  T=$(mktemp -d "${TMPDIR:-/tmp}/sp12.XXXXXX"); trap 'rm -rf "$T"' EXIT
  awk -v c="$CHR" -v s="$S" -v e="$E" '$1 == c && $3 >= s && $3 <= e' "$OUT/positions_noSNP.bed" > "$T/pos.bed"
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
  Rscript "$B/human/scripts/gt_spectrum12_eval_v6.R"
  exit 0
fi

if [ "$STAGE" = "plot" ]; then
  Rscript "$B/publication/gt_definition/scripts/GTDEF_spectrum12_plot.R"
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
