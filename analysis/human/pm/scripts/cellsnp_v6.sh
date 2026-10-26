#!/bin/bash
# cellsnp_v6.sh   (2026-10-10)
# Targeted cellSNP-lite run on the teratoma v6 NTC units (arm / sample / celltype x 5 BAM types),
# equivalent to CELLSNP_LITE + DOWNSTREAM_R (annotation strand mode, snp_off) of a Nextflow
# --callers cellsnp run, plus the call-set window states used by the v7 scoring:
#   both       results_*_v6/<unit>/04_filtered/cellsnp/<bam>/annotation/snp_off/   (published)
#   none       gt_eval_v6/prewindow/results_<tier>/...  *_step4_cov_edit_filtered.tsv.gz
#   mismatch   gt_eval_v6/prewindow/results_<tier>/...  *_step5_nomm.tsv.gz
#   clustered  results_*_v6_nomm/...                    (--mismatch_window_filter 0)
# The cellSNP-lite command is the CELLSNP_LITE process (12 threads, qual = MAPQ per BAM type).
# DOWNSTREAM_R arguments are taken from the same unit's reditools3 annotation/snp_off task
# (.command.sh), with only the caller, input and output names changed.
# Then scores the four v7 truth variants on the four window states with all callers + cellsnp
# into gt_eval_v7/variants/scoring/<state>/<variant>/<tier>/ and ranks (gt_rank_v7.R).
# cellsnp recal / recalq20 are unreported (as scAllele): cellSNP-lite has no base-quality
# filter, so BQSR qualities and the Q20 threshold do not reach it (recal == recalq20 calls).
#
#   bash human/pm/scripts/cellsnp_v6.sh submit
# Modes: manifest | call (array) | filter (array, CLASS=hi|mid|lo) | rank | submit
#$ -S /bin/bash
#$ -N cellsnp_v6
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
set -uo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
P=$PM/gt_eval_v6/cellsnp
SELF=$PM/scripts/cellsnp_v6.sh
CELLSNP=${CONDA_ENVS:-/opt/conda/envs}/cellsnp/bin/cellsnp-lite
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
SC=$PM/gt_eval_v7/variants/scoring
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
MODE=${1:-}
mkdir -p "$P/status" "$P/calls" "$P/runs"

# manifest.tsv: tier unit bam_label qual bam rd3_downstream_wd published_rd3_pre100 rss_gb class
do_manifest() {
  local out=$P/manifest.tsv
  printf "tier\tunit\tbam_label\tqual\tbam\trd3_wd\trd3_pub\trss_gb\tclass\n" > "$out"
  awk -F'\t' 'NR>1 && $3 ~ /\/annotation\/snp_off\// && $3 ~ /results_[^/]+\/(T[0-9]+_)?NTC_ctrl(__[^/]+)?\//' \
    "$PM/gt_eval_v6/reditools3_regen/manifest.tsv" | sort |
  while IFS=$'\t' read -r tier wd pub rss; do
    local inp cw bam q bl unit cls
    inp=$(readlink -f "$(ls "$wd"/*_reditools3.txt.gz)"); cw=$(dirname "$inp")
    bam=$(readlink -f "$(ls "$cw"/*.bam)")
    q=$(grep -oE '^ *-q [0-9]+' "$cw/.command.sh" | awk '{print $2}')
    bl=$(awk -F/ '{for(i=1;i<=NF;i++) if($i=="reditools3") print $(i+1)}' <<<"$pub")
    unit=${pub#"$PM"/results_*/}; unit=${unit%%/*}
    [ -s "$bam" ] && [ -s "$bam.bai" ] && [ -n "$q" ] || { echo "bad row: $tier $unit $bl bam=$bam q=$q" >&2; exit 1; }
    cls=$(awk -v g="$rss" 'BEGIN{print (g >= 64) ? "hi" : (g >= 24) ? "mid" : "lo"}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tier" "$unit" "$bl" "$q" "$bam" "$wd" "$pub" "$rss" "$cls" >> "$out"
  done || exit 1
  local c
  for c in hi mid lo; do
    awk -F'\t' -v c=$c 'NR==1 || $9==c' "$out" > "$P/manifest_$c.tsv"
  done
  wc -l "$P"/manifest*.tsv
}

row() { sed -n "$(( $2 + 1 ))p" "$1"; }

do_call() {
  local line tier unit bl q bam wd pub rss cls
  line=$(row "$P/manifest.tsv" "$SGE_TASK_ID")
  [ -n "$line" ] || { echo "empty manifest line $SGE_TASK_ID"; exit 1; }
  IFS=$'\t' read -r tier unit bl q bam wd pub rss cls <<<"$line"
  local d=$P/calls/$tier/$unit/$bl
  rm -rf "$d"; mkdir -p "$d"; cd "$d" || exit 1
  ln -sfn "$(readlink -f "$wd/hg38.analysisSet.fa")" hg38.analysisSet.fa
  ln -sfn "$(readlink -f "$wd/hg38.analysisSet.fa").fai" hg38.analysisSet.fa.fai
  grep -vE '^[[:space:]]*$' "$wd/chrlist_main_chrs.txt" > _chroms_req.txt
  samtools view -H "$bam" | awk -F'\t' '$1=="@SQ"{sub(/^SN:/,"",$2); print $2}' > _chroms_bam.txt
  awk 'NR==FNR{f[$1]=1; next} ($1 in f)' _chroms_bam.txt _chroms_req.txt \
    | awk 'NR==FNR{f[$1]=1; next} ($1 in f)' <(cut -f1 hg38.analysisSet.fa.fai) - > _chroms.txt
  echo "=== $(date) === call task=$SGE_TASK_ID $tier $unit $bl q=$q"
  echo "  bam=$bam"
  local t0=$(date +%s) rc
  $CELLSNP -s "$bam" -I "$unit" -O _csnp --chrom "$(paste -sd, _chroms.txt)" \
    --cellTAG None --UMItag None -f hg38.analysisSet.fa --minCOUNT 10 --minMAF 0 --minMAPQ "$q" \
    -p 12 --gzip > cellsnp.log 2>&1; rc=$?
  local t1=$(date +%s) out=${unit}_${bl}_cellsnp_variants.txt.gz n=NA
  if [ $rc -eq 0 ]; then
    gzip -cd _csnp/cellSNP.base.vcf.gz | awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {
        ad = dp = oth = 0; n = split($8, kv, ";")
        for (i = 1; i <= n; i++) { split(kv[i], p, "=")
            if (p[1] == "AD") ad = p[2] + 0; else if (p[1] == "DP") dp = p[2] + 0; else if (p[1] == "OTH") oth = p[2] + 0 }
        r = toupper($4); a = toupper($5)
        if (ad >= 1 && r ~ /^[ACGT]$/ && a ~ /^[ACGT]$/ && r != a) print $1, $2, r, a, dp - ad, ad, dp + oth
    }' | gzip > "$out"; rc=$?
    n=$(gzip -cd "$out" | wc -l)
    local pubd=$PM/$(cut -d/ -f1 <<<"${pub#"$PM"/}")/$unit/03_variants/cellsnp/$bl
    mkdir -p "$pubd" && cp -f "$out" "$pubd/"
    rm -rf _csnp
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" call "$tier" "$unit" "$bl" "$rc" "$((t1 - t0))" "$n" \
    > "$P/status/call__${tier}__${unit}__${bl}.tsv"
  echo "=== $(date) === done rc=$rc wall_s=$((t1 - t0)) rows=$n"
  exit $rc
}

stage() {
  local src=$1 dest=$2 e n
  rm -rf "$dest"; mkdir -p "$dest"
  for e in "$src"/*; do
    n=$(basename "$e")
    case "$n" in *_reditools3.txt.gz|*_reditools3_*) continue ;; esac
    if [ -L "$e" ]; then cp -P "$e" "$dest/"
    elif [ -d "$e" ]; then mkdir -p "$dest/$n"; cp -P "$e"/* "$dest/$n/" 2>/dev/null || true
    fi
  done
}

copy_published() {
  local srcdir=$1 dest=$2 prefix=$3 f
  mkdir -p "$dest"
  for f in "$srcdir/${prefix}"_{pre100.tsv.gz,filtered.txt.gz,summary.txt,stepwise.tsv,metrics.tsv}; do
    [ -e "$f" ] && cp -f "$f" "$dest/"
  done
}

n_ag() {
  zcat "$1" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="varType")v=i;if($i=="site")c=i};next}$v=="AG"{s[$c]=1}END{print length(s)+0}'
}

do_filter() {
  local CLASS=${CLASS:-lo} line tier unit bl q bam wd pub rss cls
  line=$(row "$P/manifest_${CLASS}.tsv" "$SGE_TASK_ID")
  [ -n "$line" ] || { echo "empty manifest line $SGE_TASK_ID CLASS=$CLASS"; exit 1; }
  IFS=$'\t' read -r tier unit bl q bam wd pub rss cls <<<"$line"
  local inp=$P/calls/$tier/$unit/$bl/${unit}_${bl}_cellsnp_variants.txt.gz
  [ -s "$inp" ] || { echo "missing calls $inp"; exit 1; }
  local base=${unit}_cellsnp_${bl}_annotation_snpoff
  local pubdir
  pubdir=$(dirname "$pub"); pubdir=${pubdir/\/04_filtered\/reditools3\//\/04_filtered\/cellsnp\/}
  local run=$P/runs/$CLASS/${SGE_TASK_ID}__$base
  echo "=== $(date) === filter CLASS=$CLASS task=$SGE_TASK_ID $tier $unit $bl rss_rd3=${rss}G"
  mk_cmd() {
    sed -e "s#--input [^ ]*#--input $(basename "$inp")#" \
        -e 's#--caller reditools3#--caller cellsnp#' \
        -e "s#_reditools3_${bl}_#_cellsnp_${bl}_#g" "$wd/.command.sh" |
      awk -v extra="$1" '{print} extra != "" && /^Rscript .*downstream_filter/{print "    " extra " \\"}'
  }
  stage "$wd" "$run"; ln -sfn "$inp" "$run/"
  mk_cmd "--save_intermediates 4" > "$run/run.sh"
  grep -q -- "--caller cellsnp" "$run/run.sh" && grep -q -- "--output_pre100 ${base}_pre100.tsv.gz" "$run/run.sh" \
    || { echo "command rewrite failed"; cat "$run/run.sh"; exit 1; }
  local t0=$(date +%s) rc
  ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?
  local t1=$(date +%s) nag=NA
  if [ $rc -eq 0 ]; then
    copy_published "$run" "$pubdir" "$base"
    cp "$run/run.log" "$pubdir/${base}_run.log"
    nag=$(n_ag "$pubdir/${base}_pre100.tsv.gz")
    local rel=${pubdir#"$PM"/results_*/}
    local pw=$PM/gt_eval_v6/prewindow/results_$tier/$rel
    mkdir -p "$pw"
    local s4=$run/${base}_filtered_step4_cov_edit_filtered.tsv.gz
    local s5=$run/${base}_filtered_step5_variable_regions_removed.tsv.gz
    if [ -s "$s4" ]; then cp -f "$s4" "$pw/${base}_step4_cov_edit_filtered.tsv.gz"; else : > "$pw/${base}_step4_cov_edit_filtered.tsv.gz"; fi
    if [ -s "$s5" ]; then cp -f "$s5" "$pw/${base}_step5_nomm.tsv.gz"; else : > "$pw/${base}_step5_nomm.tsv.gz"; fi
  else
    tail -30 "$run/run.log"
  fi
  echo "  published_rc=$rc n_AG_pre100=$nag dest=$pubdir"

  local rc_nomm=NA nag_nomm=NA
  if [ $rc -eq 0 ]; then
    local runn=${run}_nomm
    stage "$wd" "$runn"; ln -sfn "$inp" "$runn/"
    mk_cmd "--mismatch_window_filter 0" > "$runn/run.sh"
    ( cd "$runn" && bash run.sh > run.log 2>&1 ); rc_nomm=$?
    local res=${pubdir#"$PM"/}; res=${res%%/*}
    local ndest=$PM/${res}_nomm/${pubdir#"$PM"/"$res"/}
    if [ $rc_nomm -eq 0 ]; then
      copy_published "$runn" "$ndest" "$base"
      cp "$runn/run.log" "$ndest/${base}_run.log"
      nag_nomm=$(n_ag "$ndest/${base}_pre100.tsv.gz")
      rm -rf "$runn"
    else
      tail -30 "$runn/run.log"
    fi
    echo "  nomm_rc=$rc_nomm n_AG_pre100=$nag_nomm dest=$ndest"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" filter "$CLASS" "$tier" "$unit" "$bl" "$rc" "$nag" "$rc_nomm" "$nag_nomm" "$((t1 - t0))" \
    > "$P/status/filter__${tier}__${unit}__${bl}.tsv"
  [ $rc -eq 0 ] && rm -rf "$run"
  echo "=== $(date) === done rc=$rc nomm=$rc_nomm"
  [ $rc -eq 0 ] && [ "$rc_nomm" = 0 ]
}

do_rank() {
  cd "$B" || exit 1
  local n_ok n_all state v
  n_all=$(( $(wc -l < "$P/manifest.tsv") - 1 ))
  n_ok=$(cat "$P"/status/filter__*.tsv 2>/dev/null | awk -F'\t' '$6==0 && $8==0' | wc -l)
  echo "filter tasks ok: $n_ok of $n_all"
  [ "$n_ok" -eq "$n_all" ] || { echo "not all filter tasks succeeded; not ranking"; exit 1; }
  for state in none mismatch clustered both; do
    for v in v7 v7_mask v7_nomm v7_mask_nomm; do
      $RS human/scripts/gt_rank_v7.R "$SC/$state/$v" > "$SC/$state/$v/rank.log" || exit 1
    done
  done
  echo "done"
}

do_submit() {
  do_manifest || exit 1
  local n nh nm nl jc jh jm jl hold js jr
  n=$(( $(wc -l < "$P/manifest.tsv") - 1 ))
  nh=$(( $(wc -l < "$P/manifest_hi.tsv") - 1 ))
  nm=$(( $(wc -l < "$P/manifest_mid.tsv") - 1 ))
  nl=$(( $(wc -l < "$P/manifest_lo.tsv") - 1 ))
  jc=$(qsub -terse -N csnp_call -t 1-$n -tc 30 -pe threaded 12 -l mem=48G -l h_rt=72:00:00 "$SELF" call | cut -d. -f1)
  hold=$jc
  if [ $nh -gt 0 ]; then
    jh=$(qsub -terse -N csnp_filt_hi -hold_jid "$jc" -t 1-$nh -tc 20 -q himem -l mem=192G -l h_rt=24:00:00 -v CLASS=hi "$SELF" filter | cut -d. -f1); hold=$hold,$jh
  fi
  if [ $nm -gt 0 ]; then
    jm=$(qsub -terse -N csnp_filt_mid -hold_jid "$jc" -t 1-$nm -tc 30 -l mem=96G -l h_rt=24:00:00 -v CLASS=mid "$SELF" filter | cut -d. -f1); hold=$hold,$jm
  fi
  if [ $nl -gt 0 ]; then
    jl=$(qsub -terse -N csnp_filt_lo -hold_jid "$jc" -t 1-$nl -tc 40 -l mem=32G -l h_rt=24:00:00 -v CLASS=lo "$SELF" filter | cut -d. -f1); hold=$hold,$jl
  fi
  js=$(qsub -terse -N csnp_score -hold_jid "$hold" -t 1-12 -v SC_OUT="$SC" "$PM/scripts/run_gt_v7_variants_score.sh" | cut -d. -f1)
  jr=$(qsub -terse -N csnp_rank -hold_jid "$js" -l h_vmem=16G -l h_rt=2:00:00 "$SELF" rank)
  echo "call $jc ($n); filter hi=${jh:-none}($nh) mid=${jm:-none}($nm) lo=${jl:-none}($nl); score $js; rank $jr"
}

case $MODE in
  manifest) do_manifest ;;
  call)     do_call ;;
  filter)   do_filter ;;
  rank)     do_rank ;;
  submit)   do_submit ;;
  *) echo "usage: $0 manifest|call|filter|rank|submit"; exit 1 ;;
esac
