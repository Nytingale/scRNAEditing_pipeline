#!/bin/bash
# replicates_jacusa2_v8.sh   (SGE array, one task per JACUSA2 BAM type; 2026-10-15, Reviewer 3 replicates; from _v6)
# Replicate-aware JACUSA2 (native mode) on the v8 caller runs: call-1 with the four teratoma NTC BAMs as
# replicates (T1,T2,T3,T4), otherwise exactly the v8 per-teratoma JACUSA2 task (its .command.sh, BAMs and flags,
# work_gt_persample_v8). The per-replicate base counts (bases11..bases14) are summed into bases11 for
# downstream_filter.R; the score is JACUSA2's replicate-aware score.
# Downstream: the v8 arm-tier NTC_ctrl DOWNSTREAM_R task of the same BAM type and strand mode (work_gt_v8),
# rerun unchanged on that input; the step-5 table (nearby-mismatch pass = the main call state) is kept.
# Output: gt_eval_v8/replicates/results_jacusa2rep_arm/NTC_ctrl/04_filtered/jacusa2/<bam>/<strand>/snp_off/
#           NTC_ctrl_jacusa2_<bam>_<strand>_snpoff_{step5_nomm.tsv.gz,pre100.tsv.gz,stepwise.tsv}
#         gt_eval_v8/replicates/runs/jacusa2_rep/<bam>/ (raw replicate output with per-replicate counts)
#   qsub -t 1-4 replicates_jacusa2_v8.sh
#$ -S /bin/bash
#$ -N jac_rep8
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=4G
#$ -pe threaded 12
#$ -l h_rt=24:00:00
set -euo pipefail
PM=${SCRNA_DIR:-.}/human/pm
R=$PM/gt_eval_v8/replicates
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
export CONDA_PREFIX=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing
BAMS=(recal recalq20 sampleRmDup sampleUmiDup)
B=${BAMS[$((SGE_TASK_ID - 1))]}
RUN=$R/runs/jacusa2_rep/$B
rm -rf "$RUN"; mkdir -p "$RUN"

wd_of() {  # <results dir> <process> <tag> -> work dir of the COMPLETED task
  local h; h=$(awk -F'\t' -v p="$2" -v t="$3" '$4 == p && $5 == t && $7 == "COMPLETED" {print $2}' "$PM/$1/pipeline_trace.txt" | tail -1)
  [ -n "$h" ] || { echo "no COMPLETED $2 ($3) in $1" >&2; exit 1; }
  ls -d "$PM/${4}/${h}"*
}

list=""
for t in T1 T2 T3 T4; do
  jwd=$(wd_of results_gt_persample_v8 JACUSA2_CALL1 "${t}_NTC_ctrl:${B}" work_gt_persample_v8)
  bname=$(grep -E '^\s+[^ -][^ ]*\.bam\s*$' "$jwd/.command.sh" | tr -d ' ')
  bam=$(readlink -f "$jwd/$bname")
  [ -s "$bam" ] || { echo "missing BAM for $t $B ($jwd/$bname)"; exit 1; }
  [ $t = T1 ] && { B1=$bname; J1=$jwd; }
  grep -A12 "call-1" "$jwd/.command.sh" | grep -E '^\s+-[cmqP] ' | sed "s/^ *//" > "$RUN/flags_$t.txt"
  list=${list:+$list,}$bam
done
for t in T2 T3 T4; do cmp -s "$RUN/flags_T1.txt" "$RUN/flags_$t.txt" || { echo "flags differ T1 vs $t"; exit 1; }; done
echo "$B replicates: $list"; cat "$RUN/flags_T1.txt"

cd "$RUN"
for e in "$J1"/*; do [ -L "$e" ] && cp -P "$e" .; done
cp "$J1"/*.fai . 2>/dev/null || true
rm -f "$B1" "$B1.bai"
B1RE=${B1//./\\.}
sed -e "s#-r T1_NTC_ctrl_${B}_jacusa2.out#-r NTC_ctrl_${B}_jacusa2rep.out#" \
    -e "s#^\(\s*\)${B1RE}\s*\$#\1${list}#" \
    -e "s#gzip T1_NTC_ctrl_${B}_jacusa2.out#gzip NTC_ctrl_${B}_jacusa2rep.out#" \
    -e "s#scallele_io.py bam_chroms ${B1RE}#scallele_io.py bam_chroms ${list%%,*}#" \
    "$J1/.command.sh" > call.sh
grep -q "$list" call.sh || { echo "BAM list not substituted"; exit 1; }
bash call.sh > call.log 2>&1

# sum replicate counts into bases11 (downstream_filter.R jacusa2 parser)
zcat "NTC_ctrl_${B}_jacusa2rep.out.gz" | awk -F'\t' -v OFS='\t' '
  /^##/ { print; next }
  /^#/  { for (i = 1; i <= NF; i++) { if ($i ~ /^bases1[0-9]+$/) { b[++nb] = i } }
          first = b[1]; o = ""
          for (i = 1; i <= NF; i++) { skip = 0; for (k = 2; k <= nb; k++) if (i == b[k]) skip = 1
            if (!skip) o = o (o == "" ? "" : OFS) $i }
          print o; next }
  { split("0,0,0,0", s, ",")
    for (k = 1; k <= nb; k++) { split($b[k], c, ","); for (j = 1; j <= 4; j++) s[j] += c[j] }
    $first = s[1] "," s[2] "," s[3] "," s[4]
    o = ""
    for (i = 1; i <= NF; i++) { skip = 0; for (k = 2; k <= nb; k++) if (i == b[k]) skip = 1
      if (!skip) o = o (o == "" ? "" : OFS) $i }
    print o }' | gzip > "NTC_ctrl_${B}_jacusa2.out.gz"
echo "replicate calls: $(zcat "NTC_ctrl_${B}_jacusa2.out.gz" | grep -vc '^#')"

for s in annotation caller; do
  base=NTC_ctrl_jacusa2_${B}_${s}_snpoff
  dwd=$(wd_of results_gt_v8 DOWNSTREAM_R "NTC_ctrl:jacusa2:${B}:${s}:snp_off" work_gt_v8)
  grep -q -- '--mismatch_window_filter 1' "$dwd/.command.sh" || { echo "arm downstream not in the main state: $dwd"; exit 1; }
  d=$RUN/downstream_$s; rm -rf "$d"; mkdir -p "$d"
  for e in "$dwd"/*; do
    n=$(basename "$e")
    if [ -L "$e" ]; then cp -P "$e" "$d/"
    elif [ -d "$e" ]; then mkdir -p "$d/$n"; cp -P "$e"/* "$d/$n/"; fi
  done
  rm -f "$d/NTC_ctrl_${B}_jacusa2.out.gz"
  ln -s "$RUN/NTC_ctrl_${B}_jacusa2.out.gz" "$d/NTC_ctrl_${B}_jacusa2.out.gz"
  ( cd "$d" && bash "$dwd/.command.sh" > run.log 2>&1 )
  dest=$R/results_jacusa2rep_arm/NTC_ctrl/04_filtered/jacusa2/$B/$s/snp_off
  mkdir -p "$dest"; cp "$d/${base}_pre100.tsv.gz" "$d/${base}_stepwise.tsv" "$dest/"
  s5=$d/${base}_filtered_step5_variable_regions_removed.tsv.gz
  if [ -s "$s5" ]; then cp "$s5" "$dest/${base}_step5_nomm.tsv.gz"; else : > "$dest/${base}_step5_nomm.tsv.gz"; fi
  echo "$base: step5 $(zcat -f "$dest/${base}_step5_nomm.tsv.gz" | tail -n +2 | wc -l) rows, pre100 $(zcat "$dest/${base}_pre100.tsv.gz" | tail -n +2 | wc -l) rows"
done
echo done
