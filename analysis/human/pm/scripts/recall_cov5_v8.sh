#!/bin/bash
# recall_cov5_v8.sh   STAGE=manifest | submit | run (SGE array, M=<manifest subset>)   (2026-10-21)
# Caller floors for the downstream-filter sweep (user 2026-10-21: coverage 5 / 10 / 15, variant reads 2 / 3 / 4,
# editing level 1 / 2 / 3%). GATK, JACUSA2, cellSNP-lite and scAllele apply min_cov 10 inside the caller task (GATK:
# DP >= 10 when the VCF is tabulated; JACUSA2 -c; cellSNP-lite --minCOUNT; scAllele --DP), and scAllele also --AC 3,
# so their published outputs cannot give the 5-read or the 2-variant-read setting. Every caller task feeding an NTC
# task of gt_eval_v8/filter_sweep/manifest.tsv is replayed from its own .command.sh (same BAM, options and conda
# environment as the Nextflow task) with min_cov 5 (scAllele: --DP 5 --AC 2; --AB stays 0.01, the lowest swept
# editing level). The GATK unfiltered VCFs were not kept, so HaplotypeCaller is re-run. BCFtools and REDItools3 have
# no caller floor (their published outputs already hold every coverage).
#   manifest  R/manifest.tsv  id, kind (gatk | jacusa2 | cellsnp | scchrom | scgather), old_dir, pe, mem (G)
#   run       -> R/calls/<old work-dir hash path>/ (the task's outputs, same file names); R/status/<id>.tsv  id, exit
#   scgather  stages the chromosome outputs of the replayed SCALLELE_CHROM tasks (R/calls/...) instead of the old ones
#   STAGE=manifest bash recall_cov5_v8.sh; STAGE=submit bash recall_cov5_v8.sh
#$ -S /bin/bash
#$ -N recall5
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_rt=48:00:00
set -uo pipefail
export LC_ALL=C
PM=${SCRNA_DIR:-.}/human/pm
SELF=$PM/scripts/recall_cov5_v8.sh
W=$PM/gt_eval_v8/filter_sweep
R=$PM/gt_eval_v8/filter_sweep2/recall
M=${M:-$R/manifest.tsv}
STAGE=${STAGE:-run}
new_of() { echo "$R/calls/$(basename "$(dirname "$1")")/$(basename "$1")"; }

case "$STAGE" in
manifest)
  set -e; mkdir -p "$R/sub"
  awk -F'\t' 'NR>1 && $5 ~ /^(gatk|jacusa2|cellsnp|scallele)$/ {print $5"\t"$3}' "$W/manifest.tsv" |
  while IFS=$'\t' read -r c wd; do
    i=$(grep -oP -- '--input \K\S+' "$wd/.command.sh"); d=$(dirname "$(readlink "$wd/$i")")
    if [ "$c" = scallele ]; then
      for x in "$d"/*.snvs.txt; do printf "scchrom\t%s\n" "$(dirname "$(readlink "$x")")"; done
      printf "scgather\t%s\n" "$d"
    else printf "%s\t%s\n" "$c" "$d"; fi
  done | sort -u > /tmp/recall5_rows.tsv
  printf "id\tkind\told_dir\tpe\tmem\n" > "$M"
  while IFS=$'\t' read -r k d; do
    pe=$(sed -n 's/^#\$ -pe threaded \([0-9]*\).*/\1/p' "$d/.command.run")
    mem=$(sed -n 's/^#\$ -cwd -l mem=\([0-9]*\)G.*/\1/p' "$d/.command.run")
    printf "%s\t%s\t%s\t%s\n" "$k" "$d" "${pe:-1}" "${mem:-8}"
  done < /tmp/recall5_rows.tsv | awk -F'\t' -v OFS='\t' '{print NR, $0}' >> "$M"
  for c in $(awk -F'\t' 'NR>1{print $2"_"$4"_"$5}' "$M" | sort -u); do
    IFS=_ read -r k pe mem <<< "$c"
    { head -1 "$M"; awk -F'\t' -v k="$k" -v p="$pe" -v m="$mem" 'NR>1 && $2==k && $4==p && $5==m' "$M"; } > "$R/sub/${k}_${pe}_${mem}.tsv"
  done
  awk -F'\t' 'NR>1{print $2, "pe", $4, "mem", $5}' "$M" | sort | uniq -c
  ;;
submit)
  cd "$PM/../.."; chrom=""; out=""
  for f in "$R"/sub/*.tsv; do
    IFS=_ read -r k pe mem <<< "$(basename "$f" .tsv)"; [ "$k" = scgather ] && continue
    n=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N rc5_${k} -t 1-$n -tc 300 -pe threaded "$pe" -l h_rss=${mem}G,mem_free=${mem}G -l mem=${mem}G \
          -v STAGE=run,M=$f "$SELF" | cut -d. -f1)
    [ "$k" = scchrom ] && chrom=$chrom${chrom:+,}$j; out="$out $k:$j"
  done
  for f in "$R"/sub/scgather_*.tsv; do
    IFS=_ read -r k pe mem <<< "$(basename "$f" .tsv)"; n=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N rc5_gather -t 1-$n -hold_jid "$chrom" -pe threaded "$pe" -l h_rss=${mem}G,mem_free=${mem}G \
          -l mem=${mem}G -v STAGE=run,M=$f "$SELF" | cut -d. -f1); out="$out scgather:$j"
  done
  echo "$out"
  ;;
run)
  mkdir -p "$R/status"
  IFS=$'\t' read -r id kind old pe mem < <(sed -n "$(( SGE_TASK_ID + 1 ))p" "$M")
  [ "$(cut -f3 "$R/status/$id.tsv" 2>/dev/null)" = 0 ] && { echo "skip (done): $id"; exit 0; }
  new=$(new_of "$old"); rm -rf "$new"; mkdir -p "$new"
  cp "$old/.command.sh" "$new/"
  for e in "$old"/*; do
    n=$(basename "$e")
    if [ -L "$e" ]; then
      t=$(readlink "$e")
      if [ "$kind" = scgather ]; then ln -s "$(new_of "$(dirname "$t")")/$n" "$new/$n"; else cp -P "$e" "$new/"; fi
    elif [ -d "$e" ]; then mkdir -p "$new/$n"; cp -P "$e"/* "$new/$n/"
    elif [[ "$n" == *.fai || "$n" == *.dict ]]; then cp "$e" "$new/"
    fi
  done
  case $kind in
    gatk)     sed -i 's/if (cov >= 10) print/if (cov >= 5) print/' "$new/.command.sh"; pat='cov >= 5)' ;;
    jacusa2)  sed -i 's/^    -c 10 \\$/    -c 5 \\/' "$new/.command.sh"; pat='-c 5 ' ;;
    cellsnp)  sed -i 's/--minCOUNT 10 /--minCOUNT 5 /' "$new/.command.sh"; pat='--minCOUNT 5 ' ;;
    scchrom)  sed -i 's/\(scallele_chrom\.sh \S* \S* \S*\) 10 3 /\1 5 2 /' "$new/.command.sh"; pat=' 5 2 ' ;;
    scgather) pat='' ;;
  esac
  rc=0
  [ -n "$pat" ] && ! grep -qF -- "$pat" "$new/.command.sh" && { echo "floor not substituted in $new/.command.sh" >&2; rc=2; }
  if [ $rc -eq 0 ]; then
    env=$(sed -n 's/.*\/bin\/activate \(\S*\).*/\1/p' "$old/.command.run" | head -1)
    ( set +u; source ${CONDA_BASE:-/opt/conda}/bin/activate "$env"; set -u
      export PATH="$PATH:${SCRNA_DIR:-.}/nextflow/bin"
      cd "$new" && bash .command.sh > .command.out 2> .command.err ); rc=$?
  fi
  printf "%s\t%s\t%s\n" "$id" "$kind" "$rc" > "$R/status/$id.tsv"
  exit $rc
  ;;
esac
