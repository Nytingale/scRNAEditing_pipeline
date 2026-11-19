#!/bin/bash
# umi_seed_reps.sh DS=teratoma|aorta SEED=<n> STAGE=submit|split|tsetup|v9truth|v9ucov|v9collect|aortav9|nf   (user 2026-11-11)
# umi_tools seed replicates: the same umi_tools BAMs, truth and BCFtools calls as production, from BAMs deduplicated
# with --random-seed=<n> (human/pm/scripts/umi_seed_bams.sh, aorta_truth_dedup.sh), to see whether configuration
# rankings depend on the seed. Seed 8787 is the production run (in place); this script is for the other seeds.
# Only BCFtools is re-called (user): every other caller's best configuration is >= 0.02 F1 below #1, ~10x the
# deduplication differences. Only the scored units: teratoma NTC arm, T1-T4 NTC, NTC cell types; aorta 4 control units.
#   teratoma  R = human/pm/umi_seed_reps/seed<n>/
#     bams/            per-teratoma NTC / ADAR1_KO / ADARB1_KO + merged arms (umi_seed_bams.sh SEED=<n>)
#     split            NTC cell types -> celltype_bams/ (split_gt_bams_celltype.sh, plan = _celltype_frags_umimd with these BAMs)
#     truth            run_gt_truth_tso.sh scan/cands/pile/eval/truth with GT_SP12_DIR / GT_UMI_BAMDIR / GT_AGTC_DIR /
#                      GT_V7_* -> spectrum12/, agtc_test/, truth_v7*/; then gt_v9_build.sh truth/ucov/collect with
#                      GT_V9_DIR=v9/ (unit coverage on these BAMs, v8 parity check off)
#     nf               results_<tier>/ (BCFtools, --dedup umi, sampleUmiDup), step4 / step5 tables collected
#   aorta     BAMs  = mouse/aorta_gse255661/truth/umi_seed_reps/seed<n>/ (aorta_truth_dedup.sh SEED=<n> OUT=...)
#             R     = mouse/aorta_gse255661/umi_seed_reps/seed<n>/: spectrum12/ (aorta_truth_spectrum12.sh, regions and MGP
#                     SNP BED reused), eval/ (aorta_truth_eval.R), eval_v9/excluded/ (aorta_truth_v9.R), results/
#   SEED=umimd (user 2026-11-11): the same truth chain on the UMI-aware MarkDuplicates BAMs (02_gt_bams_umimd,
#   02_gt_bams_celltype_umimd, pipeline_bams_umimd; links instead of split) -> umi_seed_reps/umimd/; no caller run
#   (the umimd calls exist). Tests whether umi_tools' lead depends on sharing the truth's reads.
#   SEED=coord (user 2026-11-12): the same on the coordinate (MarkDuplicates) BAMs the coordinate configurations were
#   called on (teratoma 02_gt_bams_rmdup_v6, 02_gt_bams_celltype_rmdup_v6; aorta controls results_v8/*/02_dedup, Adar1KD
#   pipeline_bams_coord from aorta_coord_bams.sh) -> umi_seed_reps/coord/; links, no caller run.
#   Submit (login node): DS=teratoma SEED=1 HOLD=<bam merge job> bash human/pm/scripts/umi_seed_reps.sh
#$ -S /bin/bash
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=72:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
AO=$B/mouse/aorta_gse255661
SELF=$PM/scripts/umi_seed_reps.sh
DS=${DS:?DS=teratoma|aorta}; SEED=${SEED:?SEED}; STAGE=${STAGE:-submit}
[ "$SEED" != 8787 ] || { echo "seed 8787 is the in-place production rebuild" >&2; exit 1; }
TAG=seed$SEED; case $SEED in umimd|coord) TAG=$SEED ;; esac
LINKED=false; case $SEED in umimd|coord) LINKED=true ;; esac
Q=$(command -v qsub || echo /cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub)
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
if [ "$DS" = teratoma ]; then
  R=$PM/umi_seed_reps/$TAG; BAMS_DIR=$R/bams
  export GT_SP12_DIR=$R/spectrum12 GT_UMI_BAMDIR=$R/bams GT_AGTC_DIR=$R/agtc_test GT_V7_TRUTH_DIR=$R/truth_v7 \
    GT_V7_VARIANTS_DIR=$R/truth_v7_variants GT_V9_DIR=$R/v9 GT_V9_UCOV_PARITY=0
  export BAMS="$R/bams/merged_NTC_ctrl.bam $R/bams/T[1-4]_NTC_ctrl.bam $R/celltype_bams/NTC_ctrl__*.bam"
else
  R=$AO/umi_seed_reps/$TAG; BAMS_DIR=$AO/truth/umi_seed_reps/$TAG; $LINKED && BAMS_DIR=$R/bams
  export AORTA_SP12_DIR=$R/spectrum12 AORTA_BAMDIR=$BAMS_DIR AORTA_EVAL_DIR=$R/eval AORTA_V9_DIR=$R/eval_v9/excluded
fi
V="DS=$DS,SEED=$SEED"
mkdir -p "$R/logs"
LOG="-o $R/logs/ -j y"

nf_run() {  # <name> <sheet>
  local name=$1 sheet=$2 out=$R/results_$1 work=$R/work_$1 run=$R/nf_run_$1
  mkdir -p "$run"; cd "$run"
  export PS1="${PS1:-}"; export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  set +u; source ~/.bashrc 2>/dev/null || true; mamba activate scRNAEditing; set -u
  export NXF_OPTS="-Xms1g -Xmx10g"
  if [ "$DS" = teratoma ]; then
    local r=$B/human/ref_annots
    REF_ARGS=(--genome_fasta "$r/hg38.analysisSet.fa" --genome_gtf "$r/hg38.ncbiRefSeq.gtf"
      --snp_file "$B/annotations/SNP_files/dbSNP156_common_GRanges.RData"
      --strand_annots "$B/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData"
      --alu_annots "$B/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData"
      --rediportal "$r/REDIportal_hg38_v2.vcf.gz" --known_sites_vcf "$r/dbSNP156_GRCh38.vcf.gz"
      --refseq_bed "$B/review/ref_annots/hg38_RefSeq.bed" --main_chroms "$r/chrlist_main_chrs.txt")
  else
    local r=$B/mouse/ref_annots
    REF_ARGS=(--genome_fasta "$r/mm39.fa" --genome_gtf "$r/mm39.ncbiRefSeq.gtf" --snp_file "$r/GRCm39_snps_GRanges.RData"
      --strand_annots "$B/annotations/annotation_files/extended_annots_mm39.RData"
      --rediportal "$r/REDIportal_mm39.vcf" --known_sites_vcf "$r/GRCm39_snps.vcf.gz"
      --refseq_bed "$r/mm39_RefSeq.bed" --main_chroms "$r/chrlist_main_chrs.txt")
  fi
  echo "=== $(date) === seed $SEED $DS $name: BCFtools on sampleUmiDup -> $out"
  nextflow run "$B/nextflow/main.nf" -profile sge,conda --input "$sheet" --outdir "$out" "${REF_ARGS[@]}" \
    --callers bcftools --dedup umi --strand_modes "annotation,hybrid,caller" --bam_types sampleUmiDup \
    --qual 30 --min_cov 10 --bcftools_baq true --requant_editing false --percell_analyses false \
    --mismatch_window_filter true --save_intermediates 4 -work-dir "$work" -resume
  # step-4 / step-5 tables next to each pre100 (as run_pipeline_pm_gt_v8.sh STAGE=collect)
  local n_ok=0 n_bad=0 idx; idx=$(mktemp -p "$B/tmp")
  find "$work" -mindepth 3 -maxdepth 3 -name '*_pre100.tsv.gz' -type f > "$idx"
  while read -r pub; do
    local d prefix wd=""; d=$(dirname "$pub"); prefix=$(basename "$pub" _pre100.tsv.gz)
    for c in $(grep -F "/${prefix}_pre100.tsv.gz" "$idx" | xargs -r -n1 dirname); do cmp -s "$c/${prefix}_pre100.tsv.gz" "$pub" && { wd=$c; break; }; done
    [ -n "$wd" ] || { echo "NO WORK DIR: $pub" >&2; n_bad=$((n_bad + 1)); continue; }
    for st in "step4_cov_edit_filtered:step4_cov_edit_filtered" "step5_variable_regions_removed:step5_nomm"; do
      local src="$wd/${prefix}_filtered_${st%%:*}.tsv.gz" dst="$d/${prefix}_${st##*:}.tsv.gz"
      if [ -s "$src" ]; then cp "$src" "$dst"; else : > "$dst"; fi
    done
    n_ok=$((n_ok + 1))
  done < <(find "$out" -path '*/04_filtered/*' -name '*_pre100.tsv.gz')
  rm -f "$idx"
  echo "=== $(date) === $name: collected $n_ok call sets; $n_bad without a work dir"
  [ "$n_bad" -eq 0 ]
}

sheet() {  # <src sheet> <row regex> <sed on umidup_bam> -> stdout
  local src=$1 rows=$2 e=$3 col
  col=$(head -1 "$src" | tr ',' '\n' | grep -n -x umidup_bam | cut -d: -f1)
  { head -1 "$src"; tail -n +2 "$src" | grep -E "$rows" | awk -F, -v OFS=, -v c="$col" -v e="$e" \
      '{ split(e, p, "|"); n = gsub(p[2], p[3], $c); if (n != 1) { print "no swap: " $1 > "/dev/stderr"; exit 1 } print }'; }
}

case "$DS:$STAGE" in
teratoma:submit)
  H=${HOLD:?HOLD=<umi_seed_bams merge job>}
  NREG=$(wc -l < "$PM/gt_eval_v6/scan/regions.tsv"); T=$PM/scripts/run_gt_truth_tso.sh
  TV="GT_SP12_DIR=$GT_SP12_DIR,GT_UMI_BAMDIR=$GT_UMI_BAMDIR,GT_AGTC_DIR=$GT_AGTC_DIR,GT_V7_TRUTH_DIR=$GT_V7_TRUTH_DIR,GT_V7_VARIANTS_DIR=$GT_V7_VARIANTS_DIR"
  if $LINKED; then js=$($Q -terse $LOG -N s${SEED}_links -hold_jid "$H" -v $V,STAGE=links "$SELF")
  else js=$($Q -terse $LOG -N s${SEED}_split -pe threaded 8 -l mem=16G -hold_jid "$H" -v $V,STAGE=split "$SELF"); fi
  j0=$($Q -terse $LOG -N s${SEED}_tset -hold_jid "$js" -v $V,STAGE=tsetup "$SELF")
  j1=$($Q -terse $LOG -N s${SEED}_scan -t 1-"$NREG" -hold_jid "$j0" -v $TV "$T" scan | cut -d. -f1)
  j2=$($Q -terse $LOG -N s${SEED}_cands -l h_vmem=32G -hold_jid "$j1" -v $TV "$T" cands)
  j3=$($Q -terse $LOG -N s${SEED}_pile -t 1-$(( NREG * 7 )) -hold_jid "$j2" -v $TV "$T" pile | cut -d. -f1)
  j4=$($Q -terse $LOG -N s${SEED}_eval -l h_vmem=96G -hold_jid "$j3" -v $TV "$T" eval)
  j5=$($Q -terse $LOG -N s${SEED}_truth -l h_vmem=96G -hold_jid "$j4" -v $TV "$T" truth)
  j6=$($Q -terse $LOG -N s${SEED}_v9 -l h_vmem=32G -hold_jid "$j5,$js" -v $V,STAGE=v9truth "$SELF")
  j7=$($Q -terse $LOG -N s${SEED}_ucov -t 1-56 -l h_vmem=4G -hold_jid "$j6" -v $V,STAGE=v9ucov "$SELF" | cut -d. -f1)
  j8=$($Q -terse $LOG -N s${SEED}_ucol -l h_vmem=16G -hold_jid "$j7" -v $V,STAGE=v9collect "$SELF")
  jn=none; $LINKED || jn=$($Q -terse $LOG -N s${SEED}_nf -l h_vmem=64G -l h_rt=336:00:00 -hold_jid "$js" -v $V,STAGE=nf "$SELF")
  echo "seed $SEED teratoma: split $js | truth setup $j0 scan $j1 cands $j2 pile $j3 eval $j4 truth $j5 v9 $j6 ucov $j7 collect $j8 | nf $jn"
  ;;
teratoma:links)
  mkdir -p "$R/bams" "$R/celltype_bams"
  SRC=$PM/02_gt_bams_umimd; SRC_CT=$PM/02_gt_bams_celltype_umimd
  [ "$SEED" = coord ] && { SRC=$PM/02_gt_bams_rmdup_v6; SRC_CT=$PM/02_gt_bams_celltype_rmdup_v6; }
  for u in merged_NTC_ctrl merged_ADAR1_KO merged_ADARB1_KO T1_NTC_ctrl T2_NTC_ctrl T3_NTC_ctrl T4_NTC_ctrl; do
    for x in bam bam.bai; do [ -s "$SRC/$u.$x" ] || { echo "ERROR: missing $SRC/$u.$x" >&2; exit 1; }; ln -sfn "$SRC/$u.$x" "$R/bams/"; done; done
  for b in "$SRC_CT"/NTC_ctrl__*.bam; do ln -sfn "$b" "$R/celltype_bams/"; ln -sfn "$b.bai" "$R/celltype_bams/"; done
  ls "$R"/celltype_bams/*.bam | wc -l
  ;;
teratoma:split)
  P=$R/_celltype_frags; rm -rf "$P"; mkdir -p "$P" "$R/celltype_bams"
  S=$PM/_celltype_frags_umimd
  cp "$S/_barcode_prefix_map.tsv" "$S/_cells_NTC_ctrl.tsv" "$S/_keep_NTC_ctrl.txt" "$P/"
  sed "s|$PM/02_gt_bams_umimd/|$R/bams/|" "$S/_plan_NTC_ctrl.tsv" > "$P/_plan_NTC_ctrl.tsv"
  [ "$(grep -c "$R/bams/T[1-4]_NTC_ctrl.bam" "$P/_plan_NTC_ctrl.tsv")" -eq 4 ] || { echo "ERROR: plan swap" >&2; exit 1; }
  echo NTC_ctrl > "$P/_expected.txt"
  SGE_TASK_ID=1 PLANDIR="$P" OUTBAM="$R/celltype_bams" OUTBC="$P/barcodes" FRAGDIR="$P" bash $PM/scripts/split_gt_bams_celltype.sh
  [ -e "$P/_ok_NTC_ctrl" ] || { echo "ERROR: split sentinel missing" >&2; exit 1; }
  ls "$R"/celltype_bams/NTC_ctrl__*.bam | wc -l
  ;;
teratoma:tsetup)
  for u in merged_NTC_ctrl merged_ADAR1_KO merged_ADARB1_KO T1_NTC_ctrl T2_NTC_ctrl T3_NTC_ctrl T4_NTC_ctrl; do
    [ -s "$R/bams/$u.bam" ] && [ -s "$R/bams/$u.bam.bai" ] || { echo "ERROR: missing $R/bams/$u.bam" >&2; exit 1; }; done
  mkdir -p "$GT_SP12_DIR"/{scan,pileups,eval,eval_depth}
  ;;
teratoma:v9truth)  STAGE=truth bash $PM/scripts/gt_v9_build.sh ;;
teratoma:v9ucov)   STAGE=ucov bash $PM/scripts/gt_v9_build.sh ;;
teratoma:v9collect)
  STAGE=collect bash $PM/scripts/gt_v9_build.sh
  ;;
teratoma:nf)
  E="|$PM/02_gt_bams_umi_v6/|$R/bams/|"
  sheet $PM/samplesheet_pm_gt_merged_dedupv6.csv '^NTC_ctrl,' "$E" > "$R/sheet_arm.csv"
  sheet $PM/samplesheet_pm_gt_persample_dedupv6.csv '^T[1-4]_NTC_ctrl,' "$E" > "$R/sheet_sample.csv"
  { head -1 $PM/samplesheet_pm_gt_celltype_as_sample_dedupv6.csv
    for b in "$R"/celltype_bams/NTC_ctrl__*.bam; do u=$(basename "$b" .bam)
      sheet $PM/samplesheet_pm_gt_celltype_as_sample_dedupv6.csv "^${u}," "|$PM/02_gt_bams_celltype_umi_v6/|$R/celltype_bams/|" | tail -n +2; done; } > "$R/sheet_celltype.csv"
  wc -l "$R"/sheet_*.csv
  for t in arm sample celltype; do ( nf_run "$t" "$R/sheet_$t.csv" ) & done; wait
  ;;
aorta:submit)
  H=${HOLD:?HOLD=<aorta_truth_dedup merge job>}
  S12=$AO/scripts/aorta_truth_spectrum12.sh
  if $LINKED; then H=$($Q -terse $LOG -N a${SEED}_links -hold_jid "$H" -v $V,STAGE=links "$SELF"); fi
  mkdir -p "$AORTA_SP12_DIR"/{scan,pileups,snp}
  cp "$AO/truth/spectrum12/regions.tsv" "$AORTA_SP12_DIR/"
  for f in "$AO"/truth/spectrum12/snp/*; do ln -sfn "$f" "$AORTA_SP12_DIR/snp/"; done
  NREG=$(wc -l < "$AORTA_SP12_DIR/regions.tsv")
  AV="AORTA_SP12_DIR=$AORTA_SP12_DIR,AORTA_BAMDIR=$AORTA_BAMDIR,AORTA_EVAL_DIR=$AORTA_EVAL_DIR"
  j1=$($Q -terse $LOG -N a${SEED}_scan -t 1-"$NREG" -l h_rt=96:00:00 -hold_jid "$H" -v $AV "$S12" scan | cut -d. -f1)
  j2=$($Q -terse $LOG -N a${SEED}_cands -l h_vmem=32G -hold_jid "$j1" -v $AV "$S12" cands)
  j3=$($Q -terse $LOG -N a${SEED}_pile -t 1-$(( NREG * 8 )) -tc 150 -l h_rt=96:00:00 -hold_jid "$j2" -v $AV "$S12" pile | cut -d. -f1)
  j4=$($Q -terse $LOG -N a${SEED}_eval -l h_vmem=96G -hold_jid "$j3" -v $AV "$S12" eval)
  j5=$($Q -terse $LOG -N a${SEED}_v9 -l h_vmem=32G -hold_jid "$j4" -v $V,STAGE=aortav9 "$SELF")
  jn=none; $LINKED || jn=$($Q -terse $LOG -N a${SEED}_nf -l h_vmem=64G -l h_rt=336:00:00 -hold_jid "$H" -v $V,STAGE=nf "$SELF")
  echo "seed $SEED aorta: scan $j1 cands $j2 pile $j3 eval $j4 v9 $j5 | nf $jn"
  ;;
aorta:links)
  mkdir -p "$BAMS_DIR"
  for u in merged_CTRL merged_KD $(awk -F'\t' 'NR > 1 { print $1 }' "$AO/truth/design.tsv"); do
    if [ "$SEED" = umimd ]; then src=$AO/pipeline_bams_umimd/$u.bam
    elif [ "$u" = merged_CTRL ]; then src=$AO/results_v8/CTRL/02_dedup/CTRL_rmDup.bam
    elif [[ $u == Control_* ]]; then src=$AO/results_v8/$u/02_dedup/${u}_rmDup.bam
    else src=$AO/pipeline_bams_coord/$u.bam; fi
    for x in "" .bai; do [ -s "$src$x" ] || { echo "ERROR: missing $src$x" >&2; exit 1; }; ln -sfn "$src$x" "$BAMS_DIR/$u.bam$x"; done; done
  ;;
aorta:aortav9) $RS $AO/scripts/aorta_truth_v9.R ;;
aorta:nf)
  sheet $AO/pipeline_bams/samplesheet_aorta_ctrl.csv '.' "|$AO/truth/bams/|$BAMS_DIR/|" > "$R/sheet_aorta.csv"
  wc -l "$R/sheet_aorta.csv"
  nf_run aorta "$R/sheet_aorta.csv"
  ;;
*) echo "unknown DS:STAGE $DS:$STAGE" >&2; exit 1 ;;
esac
