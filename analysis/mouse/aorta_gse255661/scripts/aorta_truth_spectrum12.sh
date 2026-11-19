#!/bin/bash
# =============================================================================
# aorta_truth_spectrum12.sh  [submit|snpbed|scan|cands|pile|eval]   (eval = scripts/aorta_truth_eval.R)
#
# GSE255661 aorta ground truth, caller-independent stages. Port of
# human/pm/scripts/run_gt_spectrum12_v6.sh (the route the teratoma v8 truth is
# built from) with the same scan floor, parser and pileup flags:
#   scan   merged_CTRL (3 control mice, Myh11+ cells, UMI-dedup per mouse), per
#          40 Mb chunk: candidates per (site, strand, alt) at cov >= 5, alt >= 1,
#          >= 1% (human/scripts/mismatch_site_count.pl scan)
#   cands  merge; position BED of all candidates; MGP SNP positions FLAGGED
#          (not removed: the SNP rule is chosen at the evaluation step)
#   pile   raw per-strand counts at every candidate position, BQ / MQ 30, for
#          merged_CTRL, merged_KD and the 6 mice, per chunk
#
# Inputs (read in place): truth/bams/*.bam (aorta_truth_dedup.sh), truth/design.tsv,
#   mouse/ref_annots/{mm39.fa,chrlist_main_chrs.txt,GRCm39_snps_GRanges.RData,REDIportal_mm39.vcf}
# Outputs: truth/spectrum12/   Logs: logs/
# Do NOT run samtools on the login node; `submit` only builds regions and calls qsub.
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/mouse/aorta_gse255661/logs/
#$ -j y
#$ -pe threaded 1
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00

set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:/usr/local/bin:/usr/bin:/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export LC_ALL=C   # one sort order for every BED passed to bedtools -sorted

B="${SCRNA_DIR:-.}"
A="$B/mouse/aorta_gse255661"
SELF="$A/scripts/aorta_truth_spectrum12.sh"
OUT="${AORTA_SP12_DIR:-$A/truth/spectrum12}"   # seed replicates: AORTA_SP12_DIR / AORTA_BAMDIR (umi_seed_reps.sh)
BAMDIR="${AORTA_BAMDIR:-$A/truth/bams}"
REF="$B/mouse/ref_annots"
FASTA="$REF/mm39.fa"
MAINCHR="$REF/chrlist_main_chrs.txt"
SNP_RDATA="$REF/GRCm39_snps_GRanges.RData"
REDI="$REF/REDIportal_mm39.vcf"
REGIONS="$OUT/regions.tsv"
PARSER="$B/human/scripts/mismatch_site_count.pl"
CHUNK_BP=40000000
UNITS=(merged_CTRL merged_KD $(awk -F'\t' 'NR > 1 { print $1 }' "$A/truth/design.tsv"))

mkdir -p "$OUT"/{scan,pileups,snp}
cd "$A"
STAGE="${1:-submit}"

make_regions() {
  awk -v W="$CHUNK_BP" 'NR == FNR { keep[$1] = 1; next }
    ($1 in keep) {
      n = int(($2 + W - 1) / W); if (n < 2) n = 2
      step = int(($2 + n - 1) / n)
      for (i = 0; i < n; i++) {
        s = i * step + 1; e = (i + 1) * step; if (e > $2) e = $2
        printf "%s_%02d\t%s\t%d\t%d\t%s:%d-%d\n", $1, i + 1, $1, s, e, $1, s, e
      }
    }' "$MAINCHR" "${FASTA}.fai" > "$REGIONS"
}

if [ "$STAGE" = "submit" ]; then
  [ "${#UNITS[@]}" -eq 8 ] || { echo "ERROR: expected 8 units, got ${#UNITS[@]}" >&2; exit 1; }
  for u in "${UNITS[@]}"; do [ -s "$BAMDIR/$u.bam" ] && [ -s "$BAMDIR/$u.bam.bai" ] || { echo "ERROR: missing $u.bam(.bai)" >&2; exit 1; }; done
  make_regions; NREG=$(wc -l < "$REGIONS"); npile=$(( NREG * ${#UNITS[@]} ))
  jn=$(qsub -terse -N aor12_snp -l h_vmem=48G -l h_rt=12:00:00 "$SELF" snpbed)
  js=$(qsub -terse -N aor12_scan -t 1-"$NREG" -l h_rt=96:00:00 "$SELF" scan | cut -d. -f1)
  jc=$(qsub -terse -N aor12_cands -l h_vmem=32G -hold_jid "$jn,$js" "$SELF" cands)
  jp=$(qsub -terse -N aor12_pile -t 1-"$npile" -tc 150 -l h_rt=96:00:00 -hold_jid "$jc" "$SELF" pile | cut -d. -f1)
  echo "regions $NREG | snpbed $jn | scan $js (1-$NREG) | cands $jc | pile $jp (1-$npile; units ${UNITS[*]})"
  exit 0
fi

NREG=$(wc -l < "$REGIONS")

if [ "$STAGE" = "snpbed" ]; then
  Rscript "$B/human/scripts/dbsnp_common_to_bed.R" "$SNP_RDATA" "$OUT/snp/mgp_snps.bed.unsorted"
  sort -k1,1 -k2,2n -S 8G "$OUT/snp/mgp_snps.bed.unsorted" > "$OUT/snp/mgp_snps.bed"
  rm -f "$OUT/snp/mgp_snps.bed.unsorted"
  wc -l "$OUT/snp/mgp_snps.bed"
  exit 0
fi

if [ "$STAGE" = "scan" ]; then
  read -r CHUNK CHR S E REG <<< "$(sed -n "${SGE_TASK_ID:?}p" "$REGIONS")"
  rm -f "$OUT/scan/_ok_${CHUNK}"
  echo "=== $(date) === scan ${CHUNK} ${REG}  $(samtools --version | head -1) ==="
  samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -f "$FASTA" "$BAMDIR/merged_CTRL.bam" \
    | perl "$PARSER" scan 5 1 0.01 | gzip -c > "$OUT/scan/${CHUNK}.cand.tsv.gz"
  touch "$OUT/scan/_ok_${CHUNK}"
  echo "$CHUNK: $(zcat "$OUT/scan/${CHUNK}.cand.tsv.gz" | tail -n +2 | wc -l) candidate (site, strand, alt) rows"
  exit 0
fi

if [ "$STAGE" = "cands" ]; then
  nok=$(ls "$OUT"/scan/_ok_* 2>/dev/null | wc -l)
  [ "$nok" -eq "$NREG" ] || { echo "ERROR: $nok/$NREG scan chunks OK" >&2; exit 1; }
  [ -s "$OUT/snp/mgp_snps.bed" ] || { echo "ERROR: missing SNP BED (snpbed)" >&2; exit 1; }
  cut -f1 "$REGIONS" | while read -r c; do zcat "$OUT/scan/${c}.cand.tsv.gz" | tail -n +2; done \
    | gzip -c > "$OUT/candidates_all.tsv.gz"
  zcat "$OUT/candidates_all.tsv.gz" | awk -F'\t' -v OFS='\t' '{print $1, $2 - 1, $2}' \
    | sort -u -k1,1 -k2,2n -S 8G > "$OUT/positions_all.bed"
  bedtools intersect -u -sorted -a "$OUT/positions_all.bed" -b "$OUT/snp/mgp_snps.bed" > "$OUT/positions_mgpSNP.bed"
  bedtools intersect -v -sorted -a "$OUT/positions_all.bed" -b "$OUT/snp/mgp_snps.bed" > "$OUT/positions_noSNP.bed"
  grep -v '^#' "$REDI" | awk -F'\t' -v OFS='\t' '{print $1, $2 - 1, $2}' | sort -k1,1 -k2,2n -S 4G > "$OUT/snp/rediportal_mm39.bed"
  bedtools intersect -u -sorted -a "$OUT/positions_all.bed" -b "$OUT/snp/rediportal_mm39.bed" > "$OUT/positions_rediportal.bed"
  {
    printf "metric\tn\n"
    printf "candidate_rows\t%s\n" "$(zcat "$OUT/candidates_all.tsv.gz" | wc -l)"
    printf "positions_all\t%s\n" "$(wc -l < "$OUT/positions_all.bed")"
    printf "positions_mgpSNP\t%s\n" "$(wc -l < "$OUT/positions_mgpSNP.bed")"
    printf "positions_noSNP\t%s\n" "$(wc -l < "$OUT/positions_noSNP.bed")"
    printf "positions_rediportal\t%s\n" "$(wc -l < "$OUT/positions_rediportal.bed")"
  } | tee "$OUT/cands_stats.tsv"
  exit 0
fi

if [ "$STAGE" = "pile" ]; then
  i=$(( SGE_TASK_ID - 1 )); u=$(( i / NREG )); r=$(( i % NREG + 1 ))
  UNIT="${UNITS[$u]}"
  read -r CHUNK CHR S E REG <<< "$(sed -n "${r}p" "$REGIONS")"
  F="$OUT/pileups/${UNIT}__${CHUNK}.tsv.gz"
  rm -f "$OUT/pileups/_ok_${UNIT}__${CHUNK}"
  T="$B/tmp/aor12/${UNIT}__${CHUNK}"; rm -rf "$T"; mkdir -p "$T"; trap 'rm -rf "$T"' EXIT
  awk -v c="$CHR" -v s="$S" -v e="$E" '$1 == c && $3 >= s && $3 <= e' "$OUT/positions_all.bed" > "$T/pos.bed"
  if [ -s "$T/pos.bed" ]; then
    samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -l "$T/pos.bed" -f "$FASTA" "$BAMDIR/${UNIT}.bam" \
      | perl "$PARSER" site | gzip -c > "$F"
  else
    perl "$PARSER" site < /dev/null | gzip -c > "$F"
  fi
  touch "$OUT/pileups/_ok_${UNIT}__${CHUNK}"
  echo "$UNIT $CHUNK: $(zcat "$F" | tail -n +2 | wc -l) covered of $(wc -l < "$T/pos.bed") positions"
  exit 0
fi

if [ "$STAGE" = "eval" ]; then
  n=$(ls "$OUT"/pileups/_ok_* 2>/dev/null | wc -l); nexp=$(( NREG * ${#UNITS[@]} ))
  [ "$n" -eq "$nexp" ] || { echo "ERROR: $n/$nexp pileup chunks OK" >&2; exit 1; }
  Rscript "$A/scripts/aorta_truth_eval.R"
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
