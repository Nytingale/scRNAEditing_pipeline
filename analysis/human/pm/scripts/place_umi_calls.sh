#!/bin/bash
# place_umi_calls.sh <umi-only outdir>   (v10, user 2026-11-12)
# Called by run_pipeline_pm_gt_v8.sh / run_pipeline_aorta_v8.sh STAGE=place (UMI=1). Copies every
# <unit>/0{3_variants,4_filtered}/<caller>/sampleUmiDup subtree of a umi_tools-only run into the published tree
# (PLACE_PUB), except Cellsnp-lite hybrid stranding, which goes to the Cellsnp-lite hybrid tree (PLACE_HYB), as
# cellsnp_hybrid_v8.sh laid it out. Refuses if any destination exists (archive the old umi_tools outputs first).
set -euo pipefail
SRC=${1:?umi-only outdir}
PUB=${PLACE_PUB:?}; HYB=${PLACE_HYB:?}
[ -d "$PUB" ] || { echo "ERROR: no published tree $PUB" >&2; exit 1; }
plan=$(mktemp -p ${SCRNA_DIR:-.}/tmp); trap 'rm -f "$plan"' EXIT
find "$SRC" -mindepth 4 -maxdepth 4 -type d -name sampleUmiDup -path '*/0[34]_*/*' | sort | while read -r d; do
  rel=${d#"$SRC"/}
  if [[ "$rel" == */04_filtered/cellsnp/sampleUmiDup ]]; then
    for m in "$d"/*; do
      if [ "$(basename "$m")" = hybrid ]; then printf "%s\t%s\n" "$m" "$HYB/$rel/hybrid"
      else printf "%s\t%s\n" "$m" "$PUB/$rel/$(basename "$m")"; fi
    done
  else printf "%s\t%s\n" "$d" "$PUB/$rel"; fi
done > "$plan"
n=$(wc -l < "$plan"); [ "$n" -gt 0 ] || { echo "ERROR: nothing to place from $SRC" >&2; exit 1; }
bad=0
while IFS=$'\t' read -r s d; do [ -e "$d" ] && { echo "EXISTS: $d" >&2; bad=1; }; done < "$plan"
[ $bad -eq 0 ] || { echo "ERROR: destinations exist; nothing copied" >&2; exit 1; }
while IFS=$'\t' read -r s d; do mkdir -p "$(dirname "$d")"; cp -a "$s" "$d"; done < "$plan"
echo "placed $n subtrees from $SRC (published: $(grep -vc -F "$HYB/" "$plan"), Cellsnp-lite hybrid: $(grep -c -F "$HYB/" "$plan"))"
