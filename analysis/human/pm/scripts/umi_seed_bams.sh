#!/bin/bash
# umi_seed_bams.sh STAGE=dedup|merge   (user 2026-11-11)
# Teratoma umi_tools BAMs rebuilt with a fixed seed (--random-seed=8787), replacing 02_gt_bams_umi_v6 (unseeded:
# random_seed None in every umi_tools log, so those BAMs cannot be reproduced). Same steps as the pipeline's
# UMI_DEDUP_CHROM + UMI_DEDUP_GATHER (nextflow/modules/dedup.nf) on the TSO-clipped per-teratoma BAMs
# (samplesheet_pm_gt_persample.csv), so the output equals a pipeline run with --umi_seed 8787:
#   chunks chr1-22, X, Y (human/ref_annots/chrlist_main_chrs.txt) + _other (every other mapped contig); one umi_tools
#   run per chunk and GEM well (the seed restarts in each run), reads with UB:Z:- dropped, cell tag CR.
#   dedup  array 1-20: one task per teratoma unit -> 02_gt_bams_umi_seed8787/<T>_<ARM>.bam
#   merge  array 1-5:  one task per arm (T1-T4) -> 02_gt_bams_umi_seed8787/merged_<ARM>.bam (read count == sum checked)
#   qsub -t 1-20 -v STAGE=dedup human/pm/scripts/umi_seed_bams.sh
#   qsub -t 1-5 -v STAGE=merge -hold_jid <dedup job> human/pm/scripts/umi_seed_bams.sh
#   place  (login node, seed 8787) -> 02_gt_bams_umi_v6/ + results_gt_persample_v8/<unit>/02_dedup/ (see the stage)
#   STAGE=place bash human/pm/scripts/umi_seed_bams.sh
#$ -S /bin/bash
#$ -N umiseed
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 8
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
set -euo pipefail
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:/usr/local/bin:/usr/bin:/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:$PATH
B=${SCRNA_DIR:-.}
PM=$B/human/pm
SEED=${SEED:-8787}
# seed replicates (user 2026-11-11): SEED=<n> writes umi_seed_reps/seed<n>/bams; submit only the needed array tasks
O=$PM/02_gt_bams_umi_seed8787; [ "$SEED" = 8787 ] || O=$PM/umi_seed_reps/seed$SEED/bams
MAIN=$B/human/ref_annots/chrlist_main_chrs.txt
STAGE=${STAGE:?STAGE=dedup|merge}
NT=${NSLOTS:-8}
mkdir -p "$O"
reads() { samtools idxstats "$1" | awk '{ s += $3 + $4 } END { print s }'; }

dedup_chunk() {  # <in.bam> <chunk> <workdir>   = dedup_per_well_script with cell_tag CR
  local bam=$1 chrom=$2 w=$3/$2 out=$3/$2.dd.bam
  mkdir -p "$w"; cd "$w"
  if [ "$chrom" = "_other" ]; then
    samtools idxstats "$bam" | awk -F'\t' 'NR==FNR { m[$1] = 1; next } $1 != "*" && $3 > 0 && !($1 in m) { print $1 }' "$MAIN" - > _regions.txt
  else
    echo "$chrom" > _regions.txt
  fi
  local n_in=0 n_wells=0
  [ -s _regions.txt ] && { samtools view -b -o _in.bam "$bam" $(tr '\n' ' ' < _regions.txt); n_in=$(samtools view -c _in.bam); }
  if [ "$n_in" -eq 0 ]; then samtools view -H -b -o "$out" "$bam"; cd ..; rm -rf "$w"; echo "$chrom 0 0 0"; return; fi
  samtools view -h _in.bam | awk '/^@/ || $0 !~ /\tUB:Z:-(\t|$)/' | samtools view -b -o _clean.bam -
  mv _clean.bam _in.bam
  mkdir -p _units _dd
  samtools view -H _in.bam | awk -F'\t' '/^@RG/ { id = ""; sm = ""; for (i = 2; i <= NF; i++) { if ($i ~ /^ID:/) id = substr($i, 4); if ($i ~ /^SM:/) sm = substr($i, 4) }; print id "\t" sm }' > _rg2sm.tsv
  if [ "$(cut -f2 _rg2sm.tsv | sort -u | wc -l)" -gt 1 ]; then
    samtools split -u _units/_noRG.bam -f '_units/rg_%#.bam' _in.bam
    if [ -s _units/_noRG.bam ] && [ "$(samtools view -c _units/_noRG.bam)" -gt 0 ]; then echo "ERROR: reads without @RG ($chrom)" >&2; exit 1; fi
    rm -f _units/_noRG.bam
    for f in _units/rg_*.bam; do
      id=$(samtools view -H "$f" | awk -F'\t' '/^@RG/ { for (i = 2; i <= NF; i++) if ($i ~ /^ID:/) print substr($i, 4) }' | head -1)
      [ "$(samtools view -c "$f")" -gt 0 ] || { rm -f "$f"; continue; }
      sm=$(awk -F'\t' -v i="$id" '$1 == i { print $2 }' _rg2sm.tsv | head -1 | tr -c 'A-Za-z0-9._\n-' '_')
      mkdir -p "_units/well_$sm"; mv "$f" "_units/well_$sm/"
    done
    for d in _units/well_*; do
      x=$(basename "$d")
      if [ "$(ls "$d"/*.bam | wc -l)" -eq 1 ]; then mv "$d"/*.bam "_units/$x.bam"; else samtools merge -c -p -f -o "_units/$x.bam" "$d"/*.bam; fi
      rm -rf "$d"
    done
  else
    mv _in.bam _units/well_all.bam
  fi
  rm -f _in.bam
  for u in _units/*.bam; do
    samtools index "$u"
    umi_tools dedup --per-cell --cell-tag=CR --umi-tag=UB --extract-umi-method=tag --method=unique --random-seed=$SEED \
      -I "$u" -S "_dd/$(basename "$u")" --log="_dd/$(basename "$u" .bam).log"
  done
  n_wells=$(ls _units/*.bam | wc -l)
  if [ "$n_wells" -eq 1 ]; then mv _dd/*.bam "$out"; else samtools merge -c -p -f -o "$out" _dd/*.bam; fi
  cat _dd/*.log > "$3/$chrom.log"
  cd ..; rm -rf "$w"
  echo "$chrom $n_in $n_wells $(samtools view -c "$out")"
}
export -f dedup_chunk; export MAIN SEED

case "$STAGE" in
dedup)
  S=$(tail -n +2 "$PM/samplesheet_pm_gt_persample.csv" | sed -n "${SGE_TASK_ID}p" | cut -d, -f1)
  IN=$(tail -n +2 "$PM/samplesheet_pm_gt_persample.csv" | sed -n "${SGE_TASK_ID}p" | cut -d, -f2)
  [ -n "$S" ] && [ -s "$IN" ] || { echo "ERROR: task $SGE_TASK_ID: no unit or missing $IN" >&2; exit 1; }
  W=$B/tmp/umi_seed/s$SEED/$S; rm -rf "$W"; mkdir -p "$W"
  { cat "$MAIN"; echo _other; } | sed '/^$/d' > "$W/chunks.txt"
  echo "=== $(date) $S: $IN, $(wc -l < "$W/chunks.txt") chunks, seed $SEED ($(samtools --version | head -1); umi_tools $(umi_tools --version 2>&1 | awk '{ print $NF }'))"
  xargs -P "$NT" -I{} bash -c 'set -euo pipefail; dedup_chunk "$0" "$1" "$2"' "$IN" {} "$W" < "$W/chunks.txt" > "$W/counts.txt"
  [ "$(wc -l < "$W/counts.txt")" -eq "$(wc -l < "$W/chunks.txt")" ] || { echo "ERROR: not every chunk finished" >&2; exit 1; }
  sort -V "$W/counts.txt"
  PIECES=$(sed "s|.*|$W/&.dd.bam|" "$W/chunks.txt" | sort)
  SUM=0; for b in $PIECES; do SUM=$((SUM + $(samtools view -c "$b"))); done
  samtools merge -@ "$NT" -c -p -f -o "$O/.$S.tmp.bam" $PIECES
  samtools index -@ "$NT" "$O/.$S.tmp.bam"
  N=$(samtools view -c "$O/.$S.tmp.bam")
  [ "$N" -eq "$SUM" ] || { echo "ERROR: merged $N != sum of pieces $SUM" >&2; exit 1; }
  mv -f "$O/.$S.tmp.bam" "$O/$S.bam"; mv -f "$O/.$S.tmp.bam.bai" "$O/$S.bam.bai"
  cat $(sed "s|.*|$W/&.log|" "$W/chunks.txt") > "$O/$S.umi_dedup.log" 2>/dev/null || true
  awk -v s="$S" '{ i += $2; o += $4 } END { printf "%s\t%d\t%d\n", s, i, o }' "$W/counts.txt" > "$O/$S.counts.tsv"
  rm -rf "$W"
  echo "=== $(date) wrote $O/$S.bam ($N reads)"
  ;;
merge)
  ARMS=(AAVS1_ctrl ADAR1_KO ADARB1_KO ADARB2_KO NTC_ctrl); ARM=${ARMS[$((SGE_TASK_ID - 1))]}
  IN=(); SUM=0
  for T in T1 T2 T3 T4; do
    b=$O/${T}_$ARM.bam; [ -s "$b" ] && [ -s "$b.bai" ] || { echo "ERROR: missing $b(.bai)" >&2; exit 1; }
    IN+=("$b"); SUM=$((SUM + $(reads "$b")))
  done
  samtools merge -@ "$NT" -f -o "$O/.merged_$ARM.tmp.bam" "${IN[@]}"
  samtools index -@ "$NT" "$O/.merged_$ARM.tmp.bam"
  N=$(reads "$O/.merged_$ARM.tmp.bam")
  [ "$N" -eq "$SUM" ] || { echo "ERROR: merged $N != sum $SUM" >&2; exit 1; }
  mv -f "$O/.merged_$ARM.tmp.bam" "$O/merged_$ARM.bam"; mv -f "$O/.merged_$ARM.tmp.bam.bai" "$O/merged_$ARM.bam.bai"
  echo "=== $(date) wrote merged_$ARM.bam ($N reads = sum of T1-T4)"
  ;;
place)
  # seed 8787 into the layout the truth and caller runs read (v10 rebuild 2026-11-12): merged arms in
  # 02_gt_bams_umi_v6/, per-teratoma BAMs as results_gt_persample_v8/<unit>/02_dedup/<unit>_umiDup.bam (+ umi_tools log
  # and read counts) with links 02_gt_bams_umi_v6/<unit>.bam; refuses if any destination exists
  [ "$SEED" = 8787 ] || { echo "ERROR: place is for seed 8787" >&2; exit 1; }
  D=$PM/02_gt_bams_umi_v6
  [ ! -e "$D" ] || { echo "ERROR: $D exists" >&2; exit 1; }
  for u in $(tail -n +2 "$PM/samplesheet_pm_gt_persample.csv" | cut -d, -f1); do
    dd=$PM/results_gt_persample_v8/$u/02_dedup
    for x in bam bam.bai counts.tsv umi_dedup.log; do
      [ -e "$O/$u.$x" ] || { echo "ERROR: missing $O/$u.$x" >&2; exit 1; }
      [ ! -e "$dd/${u}_umiDup.$x" ] || { echo "ERROR: $dd/${u}_umiDup.$x exists" >&2; exit 1; }
    done
  done
  mkdir -p "$D"
  for a in AAVS1_ctrl ADAR1_KO ADARB1_KO ADARB2_KO NTC_ctrl; do mv "$O/merged_$a.bam" "$O/merged_$a.bam.bai" "$D/"; done
  for u in $(tail -n +2 "$PM/samplesheet_pm_gt_persample.csv" | cut -d, -f1); do
    dd=$PM/results_gt_persample_v8/$u/02_dedup; mkdir -p "$dd"
    for x in bam bam.bai counts.tsv umi_dedup.log; do mv "$O/$u.$x" "$dd/${u}_umiDup.$x"; done
    for x in bam bam.bai; do ln -s "$dd/${u}_umiDup.$x" "$D/$u.$x"; done
  done
  rmdir "$O"
  echo "placed: $(ls "$D"/merged_*.bam | wc -l) merged arms, $(ls "$D"/T*.bam | wc -l) per-teratoma links in $D"
  ;;
*) echo "STAGE=dedup|merge|place" >&2; exit 1 ;;
esac
