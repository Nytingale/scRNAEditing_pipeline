#!/bin/bash
# starsolo_align_runs.sh   (2026-11-18)
# STARsolo alignment of the runs in <DIR>/run_list.tsv (as written by fetch_sra_runs.sh; column 1 = run), one Slurm
# array task per run, each run aligned separately (barcodes are unique only within a GEM well). Same STAR settings as
# the pipeline's STARSOLO_ALIGN (nextflow/modules/starsolo.nf) for 10x 3' (CB_UMI_Simple, CellRanger4 TSO clipping,
# GeneFull + SJ, unique mappers, two-pass). Output: <DIR>/01_starsolobams/<run>_Aligned.sortedByCoord.out.bam (+ .bai,
# _Solo.out, logs) and <DIR>/status/align_<run>.tsv (whitelist fraction, STAR summary). Runs with a status file are skipped.
#   GENOME_PARAMS  genome params file (annotations/genomes/<genome>.yaml); its star_index is used unless STAR_INDEX is set
#   PLATFORM       10xv2 | 10xv3: CB 16 + UMI 10 / 12 and the matching whitelist in nextflow/assets/whitelists/
#   BC_READ, CDNA_READ  fastq file numbers k of fastq_files/<run>_<k>.fastq.gz (SRA numbering varies by submitter)
#   WHITELIST      overrides the platform whitelist; MIN_WL_FRACTION (default 0.5) stops on the wrong chemistry
#   The UMI length is set from the barcode read (a 26-nt 10x v3 read 1 gives a 10-nt UMI).
# Submit with SGE qsub (parallel environment "threaded"; on Slurm, an SGE-compatible qsub wrapper):
#   DIR=<dataset dir> GENOME_PARAMS=<yaml> PLATFORM=10xv3 BC_READ=2 CDNA_READ=3 NAME=<job name> TC=<max concurrent>
#     NT=<threads> MEM=<per job, e.g. 90G> [HOLD=<job id>] bash starsolo_align_runs.sh submit
#$ -S /bin/bash
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -l h_rt=72:00:00
set -euo pipefail
export LC_ALL=C
B=${SCRNA_DIR:-.}
SELF=$B/cluster/starsolo_align_runs.sh
A=${DIR:?}
mkdir -p "$A/01_starsolobams" "$A/status" "$A/logs"

case "${1:-${STAGE:-}}" in
align)
  export PATH="/usr/local/bin:/usr/bin:/bin:${PATH:-}"; export PS1="${PS1:-}"
  set +u; source ~/.bashrc 2>/dev/null || true; mamba activate scRNAEditing; set -u
  r=$(awk -F'\t' -v i="$SGE_TASK_ID" 'NR == i + 1 {print $1}' "$A/run_list.tsv")
  [ -n "$r" ] || { echo "no run for task $SGE_TASK_ID"; exit 2; }
  if [ -s "$A/status/align_$r.tsv" ]; then echo "$r already done"; exit 0; fi
  idx=${STAR_INDEX:-$(awk '$1 == "star_index:" {print $2}' "${GENOME_PARAMS:?}")}
  [ -s "$idx/SA" ] || { echo "no STAR index '$idx'"; exit 2; }
  case "${PLATFORM:?}" in
    10xv2) umi=10; min_umi=10; wl=$B/nextflow/assets/whitelists/737K-august-2016.txt ;;
    10xv3) umi=12; min_umi=10; wl=$B/nextflow/assets/whitelists/3M-february-2018.txt ;;
    *) echo "PLATFORM must be 10xv2 or 10xv3"; exit 2 ;;
  esac
  wl=${WHITELIST:-$wl}
  bc=$A/fastq_files/${r}_${BC_READ:?}.fastq.gz; cdna=$A/fastq_files/${r}_${CDNA_READ:?}.fastq.gz
  [ -s "$bc" ] && [ -s "$cdna" ] || { echo "missing $bc or $cdna"; exit 2; }
  len=$(set +o pipefail; zcat "$bc" | head -2 | tail -1 | tr -d '\n' | wc -c)
  [ "$len" -lt $((16 + umi)) ] && umi=$((len - 16))
  [ "$umi" -ge "$min_umi" ] || { echo "$r: barcode read ${len} nt is too short for $PLATFORM"; exit 2; }
  wlf=$(set +o pipefail; zcat "$bc" | head -80000 | awk 'NR % 4 == 2 {print substr($0, 1, 16)}' |
        awk -v wl="$wl" 'BEGIN { while ((getline l < wl) > 0) w[l] } { n++; if ($1 in w) k++ } END { printf "%.4f", k / n }')
  echo "$r: barcode read ${len} nt, UMI ${umi} nt, ${wlf} of sampled cell barcodes on $(basename "$wl")"
  awk -v f="$wlf" -v m="${MIN_WL_FRACTION:-0.5}" 'BEGIN { exit !(f >= m) }' ||
    { echo "$r: too few barcodes on the whitelist; wrong PLATFORM / WHITELIST?"; exit 3; }
  T=$B/tmp/starsolo_align/$r; rm -rf "$T"; mkdir -p "$T"
  cd "$A/01_starsolobams"
  echo "=== $(date) === STARsolo $r ($idx) ==="
  STAR \
    --genomeDir "$idx" \
    --readFilesIn "$cdna" "$bc" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen "$umi" \
    --soloBarcodeReadLength 0 \
    --soloCBwhitelist "$wl" \
    --clipAdapterType CellRanger4 \
    --soloFeatures GeneFull SJ \
    --twopassMode Basic \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN MD \
    --outFilterType BySJout \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM "${SORT_RAM:-20000000000}" \
    --genomeLoad NoSharedMemory \
    --outTmpDir "$T/STARtmp" \
    --outFileNamePrefix "${r}_" \
    --runThreadN "${NT:-4}"
  samtools index "${r}_Aligned.sortedByCoord.out.bam"
  um=$(awk -F'|' '/Uniquely mapped reads %/ {gsub(/[ \t]/, "", $2); print $2}' "${r}_Log.final.out")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$r" "$(date +%F_%T)" "$len" "$umi" "$wlf" "$um" > "$A/status/align_$r.tsv"
  rm -rf "$T"
  echo "=== $(date) === done $r (uniquely mapped $um) ==="
  ;;

submit)
  n=$(( $(wc -l < "$A/run_list.tsv") - 1 ))
  qsub -terse -N "${NAME:-starsolo}" -t 1-$n -tc "${TC:-1}" -pe threaded "${NT:-4}" -l mem="${MEM:-90G}" \
    ${HOLD:+-hold_jid "$HOLD"} \
    -o "$A/logs/" -e "$A/logs/" \
    -v STAGE=align,DIR="$A",GENOME_PARAMS="$(readlink -f "${GENOME_PARAMS:?}")",PLATFORM="${PLATFORM:?}",BC_READ="${BC_READ:?}",CDNA_READ="${CDNA_READ:?}",NT="${NT:-4}",WHITELIST="${WHITELIST:-}",MIN_WL_FRACTION="${MIN_WL_FRACTION:-0.5}",STAR_INDEX="${STAR_INDEX:-}",SORT_RAM="${SORT_RAM:-20000000000}" \
    "$SELF"
  ;;
*) echo "usage: DIR=... GENOME_PARAMS=... PLATFORM=... BC_READ=... CDNA_READ=... bash $0 submit"; exit 1 ;;
esac
