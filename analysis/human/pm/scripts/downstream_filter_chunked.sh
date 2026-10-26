#!/bin/bash
# downstream_filter_chunked.sh <downstream_filter*.R> <its arguments>   (2026-10-23)
# Called by the downstream filter scripts themselves for inputs >= DF_CHUNK_MIN_MB. Steps 0-6 (parse, main
# chromosomes, multi-mismatch, dbSNP / blacklist, strand, coverage / alt / edit level) are per site, so they run on
# chunks of whole chromosomes (one after another, or DF_CHUNK_WORKERS at a time); the removal of sites near other
# mismatch types (step 7) picks its window over all sites, so steps 7-9 run once on the merged step-6 tables.
# A chunk is >= DF_CHUNK_ROWS input rows (default 5e6), closed only at a chromosome change. One plain run instead
# when the input has a single chunk, a chromosome in more than one block, or a scAllele chunk exits 3.
set -euo pipefail
script=$1; shift
args=("$@")
RS=${DF_RSCRIPT:-Rscript}
ROWS=${DF_CHUNK_ROWS:-5000000}
WORKERS=${DF_CHUNK_WORKERS:-1}
export DF_CHUNK=0

get() {  # <option> -> its value ("" if absent)
  local i
  for ((i = 0; i < ${#args[@]}; i++)); do
    case "${args[i]}" in "$1") echo "${args[i + 1]:-}"; return ;; "$1"=*) echo "${args[i]#*=}"; return ;; esac
  done
}
plain() { echo "[chunked] $1: one run"; rm -rf "$tmp"; exec "$RS" "$script" "${args[@]}"; }

input=$(get --input); caller=$(get --caller); min_cov=$(get --min_cov); min_cov=${min_cov:-10}
tmp=$(mktemp -d "$PWD/.df_chunks.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

# bcftools: the same min_cov prefilter as the script's own read (header lines fail it there too)
zcat -f "$input" | awk -F'\t' -v d="$tmp" -v rows="$ROWS" -v bcf=$([ "$caller" = bcftools ] && echo 1 || echo 0) -v mc="$min_cov" '
  bcf && !($5 + 0 >= mc) { next }
  !bcf && !started && ($0 ~ /^#/ || $2 !~ /^[0-9]+$/) { hdr = hdr $0 "\n"; next }
  {
    started = 1
    if ($1 != chr) {
      if ($1 in seen) split_chr = 1
      seen[$1] = 1; chr = $1
      if (k == 0 || n >= rows) {
        if (k) close(out)
        out = sprintf("gzip -1 > %s/chunk_%04d.txt.gz", d, ++k); n = 0
        printf "%s", hdr | out
      }
    }
    print | out; n++
  }
  END { if (k) close(out); print k + 0, split_chr + 0 > (d "/split.txt") }'
read -r nk split_chr < "$tmp/split.txt"
[ "$split_chr" = 0 ] || plain "a chromosome occurs in more than one block of $input"
[ "$nk" -gt 1 ] || plain "$nk chunk(s)"
echo "[chunked] $input: $nk chunks of >= $ROWS rows, $WORKERS at a time"

# chunk arguments: chunk input, scratch outputs, no stepwise / pre100 / metrics / intermediates
drop=" --output_stepwise --output_pre100 --output_metrics --save_intermediates --input --output_filtered --output_summary "
cargs=()
for ((i = 0; i < ${#args[@]}; i++)); do
  a=${args[i]}; o=${a%%=*}
  if [[ "$drop" == *" $o "* ]]; then [[ "$a" == *=* ]] || i=$((i + 1)); continue; fi
  cargs+=("$a")
done
printf '%s\0' "${cargs[@]}" > "$tmp/cargs"

run_chunk() {  # <chunk file>
  local c=${1%.txt.gz} a=()
  while IFS= read -r -d '' x; do a+=("$x"); done < "$tmp/cargs"
  set +e
  "$RS" "$script" "${a[@]}" --input "$1" --output_filtered "$c.filtered.txt.gz" --output_summary "$c.summary.txt" \
    --save_intermediates 0 --chunk_out "$c.rds" > "$c.log" 2>&1
  echo $? > "$c.rc"
}
export -f run_chunk; export RS script tmp
ls "$tmp"/chunk_*.txt.gz | xargs -P "$WORKERS" -I{} bash -c 'run_chunk {}'

: > "$tmp/chunks.list"
for f in "$tmp"/chunk_*.txt.gz; do
  c=${f%.txt.gz}; rc=$(cat "$c.rc")
  [ "$rc" != 3 ] || plain "chunk $(basename "$c") exited 3"
  if [ "$rc" != 0 ] || [ ! -s "$c.rds" ]; then
    echo "[chunked] chunk $(basename "$c") failed (exit $rc):" >&2; tail -20 "$c.log" >&2
    [ "$rc" != 0 ] || rc=1; exit "$rc"
  fi
  echo "$c.rds" >> "$tmp/chunks.list"
done
"$RS" "$script" "${args[@]}" --chunk_in "$tmp/chunks.list"
