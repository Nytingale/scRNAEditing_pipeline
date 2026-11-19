#!/bin/bash
# fetch_sra_runs.sh   (2026-11-18)
# Download the SRA runs listed in <DIR>/run_list.tsv (tab-separated, header; column 1 = run, column 2 = ENA read_count)
# as one Slurm array task per run: prefetch -> parallel-fastq-dump --split-files --gzip (technical reads kept)
# -> <DIR>/fastq_files/<run>_<k>.fastq.gz. The spot count is checked against read_count, and the read length of each
# file goes to <DIR>/status/fetch_<run>.tsv (run, reads, date, k:len,...). Runs with a status file are skipped.
# Submit with SGE qsub (parallel environment "threaded"; on Slurm, an SGE-compatible qsub wrapper):
#   DIR=<dataset dir> NAME=<job name> TC=<max concurrent tasks> NT=<threads per task> MEM=<per job, e.g. 16G> \
#     [HOLD=<job id to wait for>] bash fetch_sra_runs.sh submit
# Keep TC low for anything that is not urgent: each task holds NT cpus for hours.
#$ -S /bin/bash
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -l h_rt=48:00:00
set -euo pipefail
export LC_ALL=C
export PATH="${SOFTWARE_DIR:-/opt/software}/sratoolkit.3.0.7-centos_linux64/bin:${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:/usr/local/bin:/usr/bin:/bin:${PATH}"
B=${SCRNA_DIR:-.}
SELF=$B/cluster/fetch_sra_runs.sh
A=${DIR:?}
F=$A/fastq_files
mkdir -p "$F" "$A/status" "$A/logs"

case "${1:-${STAGE:-}}" in
fetch)
  read -r r nexp < <(awk -F'\t' -v i="$SGE_TASK_ID" 'NR == i + 1 {print $1, $2}' "$A/run_list.tsv")
  [ -n "${r:-}" ] || { echo "no run for task $SGE_TASK_ID"; exit 2; }
  if [ -s "$A/status/fetch_$r.tsv" ]; then echo "$r already done"; exit 0; fi
  T=$B/tmp/sra_fetch/$r; mkdir -p "$T"
  echo "=== $(date) === prefetch $r ==="
  prefetch --max-size 200G -O "$T" "$r"
  sra=$(find "$T" -name "$r.sra*" -type f | head -1)
  [ -n "$sra" ] || { echo "no .sra for $r"; exit 2; }
  echo "=== $(date) === parallel-fastq-dump $r (${NT:-4} threads) ==="
  rm -f "$F/${r}"_*.fastq.gz
  parallel-fastq-dump --sra-id "$sra" --threads "${NT:-4}" --tmpdir "$T" --outdir "$F" --split-files --gzip
  for f in "$F/${r}"_*.fastq.gz; do
    [ -e "$f" ] || { echo "no fastq for $r"; exit 3; }
    k=${f##*_}; k=${k%.fastq.gz}
    n=$(( $(zcat "$f" | wc -l) / 4 ))
    [ "$n" -eq "$nexp" ] || { echo "$r file $k: $n reads vs $nexp expected"; exit 4; }
    printf "%s\t%s\n" "$k" "$(zcat "$f" | head -2 | tail -1 | tr -d '\n' | wc -c)"
  done | sort -n > "$T/readlen.tsv"
  cat "$T/readlen.tsv"
  printf "%s\t%s\t%s\t%s\n" "$r" "$nexp" "$(date +%F_%T)" "$(awk '{printf "%s%s:%s", (NR > 1 ? "," : ""), $1, $2}' "$T/readlen.tsv")" \
    > "$A/status/fetch_$r.tsv"
  rm -rf "$T"
  echo "=== $(date) === done $r ==="
  ;;

submit)
  n=$(( $(wc -l < "$A/run_list.tsv") - 1 ))
  qsub -terse -N "${NAME:-sra_fetch}" -t 1-$n -tc "${TC:-1}" -pe threaded "${NT:-4}" -l mem="${MEM:-8G}" \
    ${HOLD:+-hold_jid "$HOLD"} \
    -o "$A/logs/" -e "$A/logs/" -v STAGE=fetch,DIR="$A",NT="${NT:-4}" "$SELF"
  ;;
*) echo "usage: DIR=... bash $0 submit"; exit 1 ;;
esac
