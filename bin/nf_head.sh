#!/usr/bin/env bash
# Take the first N lines from a command (or stdin). Never fail on SIGPIPE
# (exit 141) under `set -o pipefail`.
#   nf_head.sh 10000 samtools view foo.bam
#   zcat foo.fq.gz | nf_head.sh 8
set +o pipefail
if [ $# -lt 1 ]; then
  echo "usage: nf_head.sh N [command ...]" >&2
  exit 2
fi
n=$1
shift
if [ $# -eq 0 ]; then
  head -n "$n"
else
  "$@" | head -n "$n"
fi
exit 0
