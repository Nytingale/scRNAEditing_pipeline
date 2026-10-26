#!/bin/bash
# =============================================================================
# run_gt_truth_rebuild_v6.sh
#   [init|submit|prep|scan|cands|truthpile|grid STEP|runstep STEP|sweep STEP|sweepplot STEP|freeze STEP k=v..|final]
#
# GT v6: ONE genome-wide, caller-independent truth set.
#   * site universe = every covered reference-A position (ref A on +, ref T for
#     A on -; 10x 3' forward library) in the merged NTC UMI BAM -- no callers.
#   * magnitude/depletion = merged-arm UMI BAMs; reproducibility = T1-T4 NTC
#     UMI BAMs (engine gt_eval_classify_v6.R = v3 criteria, unchanged).
#   * dbSNP156 common positions flagged, excluded by default (sweep axis).
#   * thresholds are chosen STEP BY STEP (runstep -> figures -> freeze), in the
#     order the criteria apply; each step conditioned on the steps frozen before
#     it (human/pm/gt_eval_v6/FROZEN_RULES.tsv), then applied once (`final`).
#
# Copied from run_gt_truth_rebuild_v5.sh (same PATH / samtools-1.20 / pileup
# flags -B -q30 -Q30 -d 1e6 / parser), so pileups are comparable with v5.
#
# Inputs (read in place, never copied):
#   human/pm/02_gt_bams_umi_v6/*.bam (25 = merged + T1-T4 x 5 arms; built by
#     merge_gt_umi_bams_v6.sh: per-teratoma UMI dedup, then merge)
#   human/ref_annots/hg38.analysisSet.fa(.fai), main_chrs.bed
#   human/ref_annots/REDIportal_hg38_v2.vcf.gz
#   annotations/SNP_files/dbSNP156_common_GRanges.RData
# Outputs: human/pm/gt_eval_v6/ only.   Logs: logs/
#
# Candidate floor (scan): merged NTC cov>=5, edited>=1, editing>=1% -- at or
# below the loosest value on every sweep axis, so no site that could be
# labelled at any grid point is lost. All other covered A positions are
# counted in scan/<chunk>.hist.tsv (ambiguous by construction).
#
# pe x vmem: keep h_vmem<=32G with pe=4 so jobs map to slice (not exclusive
# whole). slice has no h_rt cap -> runtimes set generously on purpose.
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
# ~/.bashrc may prepend other samtools/bedtools; re-pin the v5-parity tools.
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"

B="${SCRNA_DIR:-.}"
ENG="$B/human/scripts"
PMS="$B/human/pm/scripts"
SELF="$PMS/run_gt_truth_rebuild_v6.sh"
V6="$B/human/pm/gt_eval_v6"
FASTA="$B/human/ref_annots/hg38.analysisSet.fa"
MAINCHR="$B/human/ref_annots/main_chrs.bed"
REDI="$B/human/ref_annots/REDIportal_hg38_v2.vcf.gz"
SNP_RDATA="$B/annotations/SNP_files/dbSNP156_common_GRanges.RData"
# merged_<ARM> here = per-teratoma UMI dedup THEN merge (merge_gt_umi_bams_v6.sh);
# v5's 02_gt_bams_umi/merged_* were deduplicated after pooling T1-T4.
BAMDIR="$B/human/pm/02_gt_bams_umi_v6"
PERL="$B/misc_scripts/parse_pileup_stranded.pl"
GTDEF="$B/publication/gt_definition/scripts"

REGIONS="$V6/scan/regions.tsv"

CBED="$V6/candidates/candidates_all.bed"
NOSNP="$V6/candidates/candidates_noSNP.sites.txt"
REDI_SITES="$V6/candidates/rediportal_sites.txt"
CHUNK_BP=40000000

mkdir -p "$V6"/{scan,snp,candidates,pileups/bq30,pileups/bq20,sweeps,truth}
cd "$B"
STAGE="${1:-submit}"

# Regions: every main chromosome split into max(2, ceil(len/40Mb)) equal chunks.
make_regions() {
  awk -v W="$CHUNK_BP" 'NR==FNR { keep[$1] = 1; next }
    ($1 in keep) {
      n = int(($2 + W - 1) / W); if (n < 2) n = 2
      step = int(($2 + n - 1) / n)
      for (i = 0; i < n; i++) {
        s = i * step + 1; e = (i + 1) * step; if (e > $2) e = $2
        printf "%s_%02d\t%s\t%d\t%d\t%s:%d-%d\n", $1, i + 1, $1, s, e, $1, s, e
      }
    }' "$MAINCHR" "${FASTA}.fai" > "$REGIONS"
}

# ---- Step-wise threshold selection (2026-09-27) -----------------------------
# Thresholds are chosen in the order the criteria are applied. Each step's
# sweep is computed with every EARLIER step fixed at its frozen value
# (FROZEN_RULES.tsv) and later steps at the provisional values below, so every
# panel is conditioned on real choices. FROZEN_RULES.tsv is the single source
# of parameter values for sweep, sweepplot (GTDEF_*.R) and final.
#   FROZEN_RULES.tsv columns: step parameter value date figure
FROZEN="$V6/FROZEN_RULES.tsv"
PROV="$V6/sweeps/provisional.tsv"
STEPS=(step1_quality step2_universe step3_coverage step4_control_editing step5_reproducibility step6_depletion)

write_provisional() {
  mkdir -p "$V6/sweeps"
  printf "parameter\tvalue\nbq\t30\nsnp_filter\texcluded\nmin_cov\t10\nmin_edit_reads\t3\nmin_ctrl_edit\t0.05\nmin_cov_teratoma\t10\nconsistency_mode\tevaluable\nconsistency_min\t3\ndepletion\t0.50\nmin_cov_ko\t0\nko_min_expected\t3\n" > "$PROV"
}
step_params() {
  case "$1" in
    step1_quality)         echo "bq" ;;
    step2_universe)        echo "snp_filter" ;;
    step3_coverage)        echo "min_cov min_cov_ko ko_min_expected" ;;
    step4_control_editing) echo "min_edit_reads min_ctrl_edit" ;;
    step5_reproducibility) echo "min_cov_teratoma consistency_mode consistency_min" ;;
    step6_depletion)       echo "depletion" ;;
    *) return 1 ;;
  esac
}
is_frozen() { [ -s "$FROZEN" ] && awk -F'\t' -v k="$1" 'NR > 1 && $2 == k { f = 1 } END { exit !f }' "$FROZEN"; }
pval() {  # frozen value, else provisional
  local v=""
  is_frozen "$1" && v=$(awk -F'\t' -v k="$1" 'NR > 1 && $2 == k { v = $3 } END { print v }' "$FROZEN")
  [ -n "$v" ] || v=$(awk -F'\t' -v k="$1" '$1 == k { print $2 }' "$PROV")
  [ -n "$v" ] || { echo "ERROR: no value for $1" >&2; return 1; }
  echo "$v"
}

make_step_grid() {
  local st="$1" i p
  step_params "$st" > /dev/null || { echo "ERROR: unknown step '$st'" >&2; return 1; }
  [ -s "$PROV" ] || write_provisional
  for ((i = 0; i < ${#STEPS[@]}; i++)); do
    [ "${STEPS[$i]}" = "$st" ] && break
    for p in $(step_params "${STEPS[$i]}"); do
      is_frozen "$p" || { echo "ERROR: $st needs ${STEPS[$i]} frozen first ($p not in FROZEN_RULES.tsv)" >&2; return 1; }
    done
  done
  local BQ SNP MC MER MCE TC CM CMIN DEP
  BQ=$(pval bq); SNP=$(pval snp_filter); MC=$(pval min_cov); MER=$(pval min_edit_reads)
  MCE=$(pval min_ctrl_edit); TC=$(pval min_cov_teratoma); CM=$(pval consistency_mode)
  CMIN=$(pval consistency_min); DEP=$(pval depletion)
  local g="$V6/sweeps/$st/grid.tsv"; mkdir -p "$V6/sweeps/$st/runs"
  r() { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$@"; }
  {
    r run_id value bq snp_filter min_cov min_edit_reads min_ctrl_edit min_cov_teratoma consistency_mode consistency_min depletion
    case "$st" in
      step1_quality)  for v in 30 20; do r "bq$v" "BQ${v}_MQ${v}" "$v" "$SNP" "$MC" "$MER" "$MCE" "$TC" "$CM" "$CMIN" "$DEP"; done ;;
      step2_universe) for v in excluded included; do r "snp_$v" "$v" "$BQ" "$v" "$MC" "$MER" "$MCE" "$TC" "$CM" "$CMIN" "$DEP"; done ;;
      step3_coverage) for v in 5 6 7 8 9 10 15 20 30 50; do r "mincov_$v" "$v" "$BQ" "$SNP" "$v" "$MER" "$MCE" "$TC" "$CM" "$CMIN" "$DEP"; done ;;
      step4_control_editing)
        for k in 1 2 3 4 5; do for e in 0.01 0.02 0.05 0.10; do
          r "reads${k}_ctrl${e}" "${k};${e}" "$BQ" "$SNP" "$MC" "$k" "$e" "$TC" "$CM" "$CMIN" "$DEP"; done; done ;;
      step5_reproducibility)
        for t in 3 5 10 15 20; do for rule in strict1 strict2 strict3 evaluable; do
          if [ "$rule" = evaluable ]; then m=evaluable; n=3; else m=strict; n=${rule#strict}; fi
          r "tcov${t}_${rule}" "${t};${rule}" "$BQ" "$SNP" "$MC" "$MER" "$MCE" "$t" "$m" "$n" "$DEP"; done; done ;;
      step6_depletion)
        for d in 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90; do
          r "d$d" "$d" "$BQ" "$SNP" "$MC" "$MER" "$MCE" "$TC" "$CM" "$CMIN" "$d"; done ;;
    esac
  } > "$g"
  echo "$g"
}

# Engine call shared by `sweep` and `final`.
run_engine() {  # outdir bq snp min_cov min_edit_reads min_ctrl_edit min_cov_teratoma cmode cmin dep [min_cov_ko] [ko_min_expected]
  local out="$1" bq="$2" snp="$3" mc="$4" mer="$5" mce="$6" tc="$7" cm="$8" cmin="$9" dep="${10}" ko="${11:-$4}" kme="${12:-0}"
  local sl=()
  [ "$snp" = "excluded" ] && sl=(--site_list "$NOSNP")
  mkdir -p "$out"
  printf "bq\tsnp_filter\tmin_cov\tmin_cov_ko\tko_min_expected\tmin_edit_reads\tmin_ctrl_edit\tmin_cov_teratoma\tconsistency_mode\tconsistency_min\tdepletion\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$bq" "$snp" "$mc" "$ko" "$kme" "$mer" "$mce" "$tc" "$cm" "$cmin" "$dep" > "$out/run_params.tsv"
  Rscript "$ENG/gt_eval_classify_v6.R" \
    --truth_pileup_dir "$V6/pileups/bq${bq}" --outdir "$out" \
    --min_cov "$mc" --min_cov_ko "$ko" --ko_min_expected "$kme" --min_cov_teratoma "$tc" \
    --min_edit_reads "$mer" --min_ctrl_edit "$mce" \
    --consistency_mode "$cm" --consistency_min "$cmin" \
    --depletion "$dep" --depletion_series "$dep" \
    --rediportal_sites "$REDI_SITES" "${sl[@]}"
}

# -------------------------------------------------------------------- init ---
# Light (login-node safe): regions.tsv + sweeps/provisional.tsv.
if [ "$STAGE" = "init" ]; then
  make_regions; write_provisional
  echo "regions: $(wc -l < "$REGIONS")"; cat "$PROV"
  exit 0
fi

# -------------------------------------------------------------------- grid ---
# Light: build + print one step's grid (earlier steps must be frozen).
if [ "$STAGE" = "grid" ]; then
  g=$(make_step_grid "${2:?grid requires STEP}"); column -t "$g"
  exit 0
fi

# ----------------------------------------------------------------- runstep ---
# Login node: build a step's grid, submit its sweep array + the held figures.
if [ "$STAGE" = "runstep" ]; then
  st="${2:?runstep requires STEP}"
  g=$(make_step_grid "$st"); n=$(( $(wc -l < "$g") - 1 ))
  column -t "$g"
  rm -rf "$V6/sweeps/$st/runs"; mkdir -p "$V6/sweeps/$st/runs"
  js=$(qsub -terse -N "v6_${st%%_*}" -q slice -t 1-"$n" -pe threaded 4 -l h_vmem=32G -l h_rt=72:00:00 \
       "$SELF" sweep "$st" | cut -d. -f1)
  jp=$(qsub -terse -N "v6_${st%%_*}_plot" -q slice -l h_vmem=32G -l h_rt=48:00:00 -hold_jid "$js" \
       "$SELF" sweepplot "$st")
  echo "$st: sweep $js (array 1-$n), figures $jp"
  exit 0
fi

# ------------------------------------------------------------------ freeze ---
# Record a step's choice:  freeze STEP param=value [param=value ...]
# (FIG=<stem> overrides the figure column, default GTDEF_<STEP>.)
if [ "$STAGE" = "freeze" ]; then
  st="${2:?freeze requires STEP}"; shift 2
  step_params "$st" > /dev/null || { echo "ERROR: unknown step '$st'" >&2; exit 1; }
  [ -s "$FROZEN" ] || printf "step\tparameter\tvalue\tdate\tfigure\n" > "$FROZEN"
  for kv in "$@"; do
    k="${kv%%=*}"; v="${kv#*=}"
    echo " $(step_params "$st") " | grep -q " $k " || { echo "ERROR: $k is not a $st parameter" >&2; exit 1; }
    awk -F'\t' -v OFS='\t' -v k="$k" 'NR == 1 || $2 != k' "$FROZEN" > "$FROZEN.tmp" && mv "$FROZEN.tmp" "$FROZEN"
    printf "%s\t%s\t%s\t%s\t%s\n" "$st" "$k" "$v" "$(date +%F)" "${FIG:-GTDEF_${st}}" >> "$FROZEN"
  done
  column -t "$FROZEN"
  exit 0
fi

# ------------------------------------------------------------------ submit ---
if [ "$STAGE" = "submit" ]; then
  n=$(cd "$BAMDIR" && ls *.bam 2>/dev/null | wc -l)
  [ "$n" -eq 25 ] || { echo "ERROR: expected 25 BAMs in $BAMDIR, got $n" >&2; exit 1; }
  make_regions; write_provisional
  nreg=$(wc -l < "$REGIONS")
  echo "regions: $nreg   UMI BAMs: $n"
  j1=$(qsub -terse -N v6_prep -q slice -l h_vmem=16G -l h_rt=48:00:00 "$SELF" prep | cut -d. -f1)
  echo "1 prep      -> $j1"
  j2=$(qsub -terse -N v6_scan -q slice -t 1-"$nreg" -l h_vmem=8G -l h_rt=96:00:00 \
       "$SELF" scan | cut -d. -f1)
  echo "2 scan      -> $j2 (array 1-$nreg)"
  j3=$(qsub -terse -N v6_cands -q slice -l h_vmem=32G -l h_rt=48:00:00 \
       -hold_jid "$j1,$j2" "$SELF" cands | cut -d. -f1)
  echo "3 cands     -> $j3"
  j4=$(qsub -terse -N v6_pile -q slice -t 1-50 -l h_vmem=16G -l h_rt=72:00:00 \
       -hold_jid "$j3" "$SELF" truthpile | cut -d. -f1)
  echo "4 truthpile -> $j4 (array 1-50)"
  echo "Threshold steps are run one at a time after this: runstep step1_quality, freeze, runstep step2_universe, ..."
  exit 0
fi

# ------------------------------------------------------------------ robust ---
# Robustness lattice (7 axes, 28,125 combos; KO gate = informative-KO m; per-
# teratoma NTC coverage swept 8-12) around the fixed centre, via the
# vectorised evaluator human/scripts/gt_robustness_grid_v6.R (engine rules
# reproduced exactly; checked against real engine runs by `robustparity_*`).
# Login node: `robust` builds the grid and submits everything, chained.
ROB="$ENG/gt_robustness_grid_v6.R"
RDIR="$V6/robustness"
if [ "$STAGE" = "robust" ]; then
  mkdir -p "$RDIR"
  Rscript "$ROB" grid
  nb=$(( $(cut -f"$(head -1 "$RDIR/grid.tsv" | tr '\t' '\n' | grep -nx block | cut -d: -f1)" "$RDIR/grid.tsv" | tail -n +2 | sort -u | wc -l) ))
  np=$(( $(wc -l < "$RDIR/parity_combos.tsv") - 1 ))
  jt=$(qsub -terse -N v6_rob_task -q slice -t 1-"$nb" -l h_vmem=24G -l h_rt=24:00:00 "$SELF" robusttask | cut -d. -f1)
  jg=$(qsub -terse -N v6_rob_gather -q slice -l h_vmem=48G -l h_rt=48:00:00 -hold_jid "$jt" "$SELF" robustgather)
  jp=$(qsub -terse -N v6_rob_parity -q slice -t 1-"$np" -pe threaded 4 -l h_vmem=32G -l h_rt=72:00:00 \
       "$SELF" robustparity_run | cut -d. -f1)
  jc=$(qsub -terse -N v6_rob_pcheck -q slice -l h_vmem=24G -l h_rt=12:00:00 -hold_jid "$jt,$jp" "$SELF" robustparity_check)
  jm=$(qsub -terse -N v6_rob_mcurve -q slice -l h_vmem=24G -l h_rt=24:00:00 "$SELF" robustmcurve)
  jf=$(qsub -terse -N v6_rob_plot -q slice -l h_vmem=32G -l h_rt=48:00:00 -hold_jid "$jg,$jc,$jm" "$SELF" sweepplot robustness)
  echo "robust: tasks $jt (1-$nb)  gather $jg  parity runs $jp (1-$np)  parity check $jc  m curve $jm  figures $jf"
  exit 0
fi
if [ "$STAGE" = "robusttask" ]; then Rscript "$ROB" task "$SGE_TASK_ID"; exit 0; fi
if [ "$STAGE" = "robustgather" ]; then Rscript "$ROB" gather; exit 0; fi
if [ "$STAGE" = "robustmcurve" ]; then Rscript "$ROB" mcurve; exit 0; fi
if [ "$STAGE" = "robustparity_check" ]; then Rscript "$ROB" parity; exit 0; fi
if [ "$STAGE" = "robustparity_run" ]; then
  row=$(awk -F'\t' -v r=$(( SGE_TASK_ID + 1 )) 'NR == r' "$RDIR/parity_combos.tsv")
  IFS=$'\t' read -r CB BLK MC KO KME MER MCE DEP TC CM CMIN <<< "$row"
  echo "=== $(date) === parity engine run $CB: min_cov=$MC min_cov_ko=$KO ko_min_expected=$KME reads=$MER ctrl_edit=$MCE d=$DEP tcov=$TC $CM $CMIN ==="
  run_engine "$RDIR/parity/$CB" 30 excluded "$MC" "$MER" "$MCE" "$TC" "$CM" "$CMIN" "$DEP" "$KO" "$KME"
  exit 0
fi

# ------------------------------------------------------------------- score ---
# Score every method configuration's NTC calls against the frozen truth
# (human/scripts/gt_score_v6.R). Login node: `score TIER` submits one job;
# `scoretask TIER` is what that job runs.  TIER = arm | sample | celltype
if [ "$STAGE" = "score" ]; then
  t="${2:?score requires TIER (arm|sample|celltype)}"
  qsub -terse -N "v6_score_${t}" -q slice -l h_vmem=32G -l h_rt=24:00:00 "$SELF" scoretask "$t"
  exit 0
fi
if [ "$STAGE" = "scoretask" ]; then
  [ -s "$V6/truth/gt_truth_sites.tsv.gz" ] || { echo "ERROR: run final first" >&2; exit 1; }
  Rscript "$ENG/gt_score_v6.R" "${2:?scoretask requires TIER}"
  exit 0
fi

# ---------------------------------------------------------------- spectrum ---
# Step 1 evidence that uses no ground-truth threshold: strand-aware substitution
# spectrum of every aligned base in merged NTC and merged ADAR1_KO, at BQ/MQ 30
# and 20, over the same genome chunks as `scan`. Read-string parsing follows
# misc_scripts/parse_pileup_stranded.pl (forward), extended to C/G reference
# positions and summed per chunk (human/scripts/mismatch_spectrum_count.pl).
# Login node: `spectrum` submits the array, the gather and the step-1 figure.
SPEC="$V6/spectrum"
MMC="$ENG/mismatch_spectrum_count.pl"
if [ "$STAGE" = "spectrum" ]; then
  mkdir -p "$SPEC/chunks"; nreg=$(wc -l < "$REGIONS")
  js=$(qsub -terse -N v6_spec -q slice -t 1-"$nreg" -l h_vmem=8G -l h_rt=96:00:00 "$SELF" spectrumtask | cut -d. -f1)
  jg=$(qsub -terse -N v6_spec_gather -q slice -l h_vmem=8G -l h_rt=12:00:00 -hold_jid "$js" "$SELF" spectrumgather)
  jf=$(qsub -terse -N v6_step1_plot -q slice -l h_vmem=32G -l h_rt=48:00:00 -hold_jid "$jg" "$SELF" sweepplot step1_quality)
  echo "spectrum: chunks $js (1-$nreg)  gather $jg  figures $jf"
  exit 0
fi
if [ "$STAGE" = "spectrumtask" ]; then
  line=$(sed -n "${SGE_TASK_ID}p" "$REGIONS")
  [ -n "$line" ] || { echo "ERROR: no region for task ${SGE_TASK_ID}" >&2; exit 1; }
  read -r CHUNK CHR S E REG <<< "$line"
  out="$SPEC/chunks/${CHUNK}.tsv"; rm -f "$SPEC/chunks/_ok_${CHUNK}" "$out"
  echo "=== $(date) === v6 spectrum task ${SGE_TASK_ID}: ${CHUNK} ${REG}  $(samtools --version | head -1) ==="
  for bam in merged_NTC_ctrl merged_ADAR1_KO; do for q in 30 20; do
    samtools mpileup -B -q "$q" -Q "$q" -d 1000000 -r "$REG" -f "$FASTA" "$BAMDIR/${bam}.bam" \
      | perl "$MMC" "${bam}:${q}" >> "$out"
  done; done
  [ -s "$out" ] || { echo "ERROR: empty spectrum output for $CHUNK" >&2; exit 1; }
  touch "$SPEC/chunks/_ok_${CHUNK}"
  echo "=== $(date) === done ${CHUNK} ==="
  exit 0
fi
if [ "$STAGE" = "spectrumgather" ]; then
  nreg=$(wc -l < "$REGIONS"); nok=$(ls "$SPEC"/chunks/_ok_* 2>/dev/null | wc -l)
  [ "$nok" -eq "$nreg" ] || { echo "ERROR: $nok/$nreg spectrum chunks OK" >&2; exit 1; }
  cat "$SPEC"/chunks/*.tsv | awk -F'\t' -v OFS='\t' '
    { split($1, a, ":"); k = a[1] OFS a[2] OFS $2 OFS $3; n[k] += $4 }
    END { print "sample", "bq_mq", "rna_ref", "rna_base", "n"; for (k in n) printf "%s\t%.0f\n", k, n[k] }' \
    | { read -r h; echo "$h"; sort; } > "$SPEC/mismatch_spectrum.tsv"
  echo "spectrum: $(( $(wc -l < "$SPEC/mismatch_spectrum.tsv") - 1 )) rows -> $SPEC/mismatch_spectrum.tsv"
  exit 0
fi

# -------------------------------------------------------------------- prep ---
if [ "$STAGE" = "prep" ]; then
  echo "=== $(date) === v6 prep: dbSNP common -> BED ==="
  [ -s "$REGIONS" ] || make_regions
  [ -s "$PROV" ]    || write_provisional
  Rscript "$ENG/dbsnp_common_to_bed.R" "$SNP_RDATA" "$V6/snp/dbSNP156_common.bed"
  wc -l "$V6/snp/dbSNP156_common.bed"
  exit 0
fi

# -------------------------------------------------------------------- scan ---
if [ "$STAGE" = "scan" ]; then
  line=$(sed -n "${SGE_TASK_ID}p" "$REGIONS")
  [ -n "$line" ] || { echo "ERROR: no region for task ${SGE_TASK_ID}" >&2; exit 1; }
  read -r CHUNK CHR S E REG <<< "$line"
  BAM="$BAMDIR/merged_NTC_ctrl.bam"
  CAND="$V6/scan/${CHUNK}.cand.tsv.gz"; HIST="$V6/scan/${CHUNK}.hist.tsv"
  rm -f "$V6/scan/_ok_${CHUNK}"
  echo "=== $(date) === v6 scan task ${SGE_TASK_ID}: ${CHUNK} ${REG}  $(samtools --version | head -1) ==="
  samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -f "$FASTA" "$BAM" \
    | perl "$PERL" forward \
    | awk -F'\t' -v OFS='\t' -v cand="gzip -c > $CAND" -v hist="$HIST" '
      NR == 1 { print | cand; next }
      {
        cov = $4 + 0; alt = $5 + 0
        if (cov < 1) next
        ref = toupper($3); ncov[ref]++
        if (cov >= 5 && alt >= 1 && $16 + 0 >= 0.01) { print | cand; ncand++ }
        cb = (cov >= 100) ? 100 : (cov >= 50) ? 50 : (cov >= 30) ? 30 : (cov >= 20) ? 20 : \
             (cov >= 15) ? 15 : (cov >= 10) ? 10 : (cov >= 5) ? 5 : (cov >= 3) ? 3 : cov
        ab = (alt >= 10) ? 10 : alt
        e = $16 + 0
        lb = (e >= 0.50) ? 0.50 : (e >= 0.20) ? 0.20 : (e >= 0.10) ? 0.10 : \
             (e >= 0.05) ? 0.05 : (e >= 0.02) ? 0.02 : (e >= 0.01) ? 0.01 : 0
        h[ref OFS cb OFS ab OFS lb]++
      }
      END {
        close(cand)
        print "ref", "cov_bin_min", "edited_reads_capped10", "edit_bin_min", "n_positions" > hist
        for (k in h) print k, h[k] > hist
        printf "covered A(+)=%d  T(-)=%d  candidates=%d\n", ncov["A"], ncov["T"], ncand > "/dev/stderr"
      }'
  [ -s "$CAND" ] && [ -s "$HIST" ] || { echo "ERROR: empty scan output for $CHUNK" >&2; exit 1; }
  touch "$V6/scan/_ok_${CHUNK}"
  echo "=== $(date) === done ${CHUNK} ==="
  exit 0
fi

# ------------------------------------------------------------------- cands ---
if [ "$STAGE" = "cands" ]; then
  nreg=$(wc -l < "$REGIONS"); nok=$(ls "$V6"/scan/_ok_* 2>/dev/null | wc -l)
  [ "$nok" -eq "$nreg" ] || { echo "ERROR: $nok/$nreg scan chunks OK" >&2; exit 1; }
  [ -s "$V6/snp/dbSNP156_common.bed" ] || { echo "ERROR: missing dbSNP BED (prep)" >&2; exit 1; }
  echo "=== $(date) === v6 cands: merge chunks, flag dbSNP common + REDIportal ($(bedtools --version)) ==="
  cut -f1 "$REGIONS" | while read -r c; do
    zcat "$V6/scan/${c}.cand.tsv.gz" | tail -n +2
  done | awk -F'\t' -v OFS='\t' '{ s = (toupper($3) == "A") ? "+" : "-"; print $1, $2 - 1, $2, $1 ":" $2, ".", s }' \
    | sort -k1,1 -k2,2n -S 8G > "$CBED"
  bedtools intersect -u -a "$CBED" -b "$V6/snp/dbSNP156_common.bed" | cut -f4 > "$V6/candidates/dbsnp_common_sites.txt"
  bedtools intersect -u -a "$CBED" -b "$REDI" | cut -f4 > "$REDI_SITES"
  awk -F'\t' -v OFS='\t' -v nosnp="$NOSNP" '
    FILENAME == ARGV[1] { snp[$1] = 1; next }
    FILENAME == ARGV[2] { redi[$1] = 1; next }
    FNR == 1 { print "site", "strand", "in_dbsnp_common", "in_rediportal" }
    {
      print $4, $6, ($4 in snp) ? 1 : 0, ($4 in redi) ? 1 : 0
      if (!($4 in snp)) print $4 > nosnp
    }' "$V6/candidates/dbsnp_common_sites.txt" "$REDI_SITES" "$CBED" > "$V6/candidates/candidates_annot.tsv"
  {
    printf "metric\tn\n"
    printf "candidates_all\t%s\n"        "$(wc -l < "$CBED")"
    printf "in_dbsnp_common\t%s\n"       "$(wc -l < "$V6/candidates/dbsnp_common_sites.txt")"
    printf "candidates_noSNP\t%s\n"      "$(wc -l < "$NOSNP")"
    printf "in_rediportal\t%s\n"         "$(wc -l < "$REDI_SITES")"
  } | tee "$V6/candidates/cands_stats.tsv"
  exit 0
fi

# --------------------------------------------------------------- truthpile ---
if [ "$STAGE" = "truthpile" ]; then
  mapfile -t UNITS < <(cd "$BAMDIR" && ls *.bam | sed 's/\.bam$//' | sort)
  idx=$(( (SGE_TASK_ID - 1) % 25 )); BQ=$(( SGE_TASK_ID <= 25 ? 30 : 20 ))
  UNIT="${UNITS[$idx]}"
  [ -s "$CBED" ] || { echo "ERROR: missing $CBED" >&2; exit 1; }
  echo "=== task ${SGE_TASK_ID}: ${UNIT} (UMI) BQ${BQ}/MQ${BQ} over $(wc -l < "$CBED") candidate sites ==="
  bash "$ENG/gt_site_pileup.sh" --bam "$BAMDIR/${UNIT}.bam" --bed "$CBED" --fasta "$FASTA" \
    --out "$V6/pileups/bq${BQ}/pileup_${UNIT}.tsv" --libtype forward --min-bq "$BQ" --min-mq "$BQ"
  exit 0
fi

# ------------------------------------------------------------------- sweep ---
if [ "$STAGE" = "sweep" ]; then
  st="${2:?sweep requires STEP}"; g="$V6/sweeps/$st/grid.tsv"
  row=$(awk -F'\t' -v r=$(( SGE_TASK_ID + 1 )) 'NR == r' "$g")
  [ -n "$row" ] || { echo "ERROR: no grid row ${SGE_TASK_ID} in $g" >&2; exit 1; }
  IFS=$'\t' read -r RUN VAL BQ SNP MC MER MCE TC CM CMIN DEP <<< "$row"
  for u in merged_NTC_ctrl merged_ADAR1_KO merged_ADARB1_KO T1_NTC_ctrl T2_NTC_ctrl T3_NTC_ctrl T4_NTC_ctrl; do
    [ -s "$V6/pileups/bq${BQ}/pileup_${u}.tsv" ] || { echo "ERROR: missing pileup bq${BQ}/${u}" >&2; exit 1; }
  done
  echo "=== $(date) === v6 $st ${SGE_TASK_ID}: ${RUN} (${VAL}) ==="
  run_engine "$V6/sweeps/$st/runs/$RUN" "$BQ" "$SNP" "$MC" "$MER" "$MCE" "$TC" "$CM" "$CMIN" "$DEP"
  exit 0
fi

# --------------------------------------------------------------- sweepplot ---
if [ "$STAGE" = "sweepplot" ]; then
  st="${2:?sweepplot requires STEP (or all)}"
  echo "=== $(date) === v6 sweepplot (GTDEF) $st ==="
  Rscript "$GTDEF/GTDEF_prep_source_data.R" pm_gt "$st"
  Rscript "$GTDEF/GTDEF_gtdefinition_plot.R" pm_gt "$st"
  exit 0
fi

# ------------------------------------------------------------------- final ---
# Parameters come ONLY from FROZEN_RULES.tsv; every step must be frozen.
if [ "$STAGE" = "final" ]; then
  for st in "${STEPS[@]}"; do for p in $(step_params "$st"); do
    is_frozen "$p" || { echo "ERROR: $p ($st) not in $FROZEN -- freeze every step first" >&2; exit 1; }
  done; done
  echo "=== $(date) === v6 FINAL truth set with frozen rules ==="
  column -t "$FROZEN"
  run_engine "$V6/truth" "$(pval bq)" "$(pval snp_filter)" "$(pval min_cov)" "$(pval min_edit_reads)" \
    "$(pval min_ctrl_edit)" "$(pval min_cov_teratoma)" "$(pval consistency_mode)" "$(pval consistency_min)" "$(pval depletion)" \
    "$(pval min_cov_ko)" "$(pval ko_min_expected)"
  V6="$V6" Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))
v6 = Sys.getenv("V6")
t = fread(file.path(v6, "truth", "gt_truth_sites.tsv.gz"))
tr = t[class == "TRUE", .(site)]
tr[, c("chr", "pos") := tstrsplit(site, ":", fixed = TRUE)]
ann = fread(file.path(v6, "candidates", "candidates_annot.tsv"), select = c("site", "strand"))
tr = merge(tr, ann, by = "site", all.x = TRUE)
tr[, pos := as.integer(pos)]
setorder(tr, chr, pos)
fwrite(tr[, .(chr, start = pos - 1L, end = pos, site, score = ".", strand)],
       file.path(v6, "truth", "v6_truth_TRUE.bed"), sep = "\t", col.names = FALSE)
h = rbindlist(lapply(list.files(file.path(v6, "scan"), "\\.hist\\.tsv$", full.names = TRUE), fread))
cs = fread(file.path(v6, "candidates", "cands_stats.tsv"))
cls = t[, .N, by = class]
n_lab = sum(cls[class != "ambiguous", N])
s = rbind(
  data.table(metric = "covered_A_positions_cov_ge1", n = sum(h$n_positions)),
  data.table(metric = "covered_A_positions_cov_ge10", n = sum(h[cov_bin_min >= 10, n_positions])),
  data.table(metric = cs$metric, n = cs$n),
  data.table(metric = paste0("class_", cls$class), n = cls$N),
  # everything covered but not TRUE/FP: ambiguous candidates + dbSNP-excluded + below-floor positions
  data.table(metric = "unlabelled_covered_A_positions", n = sum(h$n_positions) - n_lab))
fwrite(s, file.path(v6, "truth", "v6_summary.tsv"), sep = "\t")
print(s)
RS
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2
exit 1
