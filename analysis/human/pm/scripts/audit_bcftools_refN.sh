#!/bin/bash
# audit_bcftools_refN.sh  [list|<array task>|collect]   (2026-10-14)
# BCFTOOLS_MPILEUP bug check: chromosomes whose variant rows all have REF = N (reference not read, most
# likely a race on building hg38.analysisSet.fa.fai in the task dir by the parallel per-chromosome
# workers). downstream_filter.R then builds varType "N<alt>", so every site on such a chromosome is lost.
#   list     -> logs/../human/pm/audit_bcftools_refN/files.tsv (all published *_bcftools_variants.txt.gz)
#   tasks    qsub -t 1-<n chunks> audit_bcftools_refN.sh   (25 files per task)
#   collect  -> audit_bcftools_refN/{per_chrom,per_file}.tsv
# per_chrom: file, chr, rows, refN rows. A chromosome is "broken" if >= 99% of its rows have REF = N.
#$ -S /bin/bash
#$ -N bcf_refN
#$ -q slice
#$ -cwd
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
set -euo pipefail
PM=${SCRNA_DIR:-.}/human/pm
# AUDIT_OUT / AUDIT_ROOTS (space-separated result dirs) override the teratoma defaults (-v for qsub).
A=${AUDIT_OUT:-$PM/audit_bcftools_refN}
ROOTS=${AUDIT_ROOTS:-"$PM/results_gt_v6 $PM/results_gt_persample_v6 $PM/results_gt_celltype_v6 $PM/results_gt_split_B"}
CH=25
mkdir -p "$A/parts"

if [ "${1:-}" = "list" ]; then
  for d in $ROOTS; do
    find -L "$d" -maxdepth 5 -path '*03_variants/bcftools*' -name '*_bcftools_variants.txt.gz'
  done | sort > "$A/files.tsv"
  n=$(wc -l < "$A/files.tsv"); echo "$n files -> $(( (n + CH - 1) / CH )) tasks"
  exit 0
fi

if [ "${1:-}" = "collect" ]; then
  { printf 'file\tchr\trows\trefN\n'; cat "$A"/parts/*.tsv; } > "$A/per_chrom.tsv"
  awk -F'\t' -v OFS='\t' 'NR > 1 {f[$1]; r[$1] += $3; if ($3 > 0 && $4 >= 0.99 * $3) {b[$1]++; bl[$1] = bl[$1] (bl[$1] ? "," : "") $2}}
    END {print "file", "rows", "n_broken_chr", "broken_chr"; for (k in f) print k, r[k], b[k] + 0, (bl[k] ? bl[k] : "-")}' \
    "$A/per_chrom.tsv" | sort > "$A/per_file.tsv"
  echo "files: $(( $(wc -l < "$A/per_file.tsv") - 1 )); with a broken chromosome: $(awk -F'\t' 'NR > 1 && $3 > 0' "$A/per_file.tsv" | wc -l)"
  exit 0
fi

t=${SGE_TASK_ID:?}
o=$A/parts/$(printf '%04d' "$t").tsv
: > "$o.tmp"
sed -n "$(( (t - 1) * CH + 1 )),$(( t * CH ))p" "$A/files.tsv" | while read -r f; do
  zcat "$f" | awk -F'\t' -v OFS='\t' -v f="$f" '{n[$1]++; if ($3 == "N") z[$1]++} END {for (c in n) print f, c, n[c], z[c] + 0}' >> "$o.tmp"
done
mv "$o.tmp" "$o"
