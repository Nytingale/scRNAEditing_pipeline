#!/bin/bash
# ============================================================================
# run_BENCHMARK_figures.sh — rebuild every figure in the benchmark-grid
# deliverable from the prepared source data.
# ============================================================================
# Assumes the source-data chain has already been run:
#   1. the Nextflow pipeline          -> summary/collated/all_results_summary.tsv
#   2. Rscript scripts/BENCHMARK_prep_source_data.R   (per-run -> aggregated)
# This script is step 3 and is cheap (seconds per panel), so it is safe to run
# interactively on a login node.
#
# >>> PLOTTING IS DELIBERATELY OUTSIDE THE PIPELINE. <<<
# COLLATE_RESULTS used to system2() the plotting script, so a ggplot bug failed
# the whole run: a 35-second OOM discarded ~28 h of compute and cost two full
# relaunches of nf_pm_ps even though every result was already computed. Figures
# are now rebuilt from the prepared tables alone, so a broken figure costs a
# re-run of THIS script and nothing else.
#
# BENCHMARK_grid_source_data.R writes the plotted tables to source_data/ and
# BENCHMARK_grid_plot.R draws each lettered panel and the combined figure they
# assemble into from those tables only. Every SVG goes to figures/; every TSV
# (source data, .sources.tsv sidecars and figure_source_map.tsv) goes to
# source_data/.
#
# Usage:
#   bash .../run_BENCHMARK_figures.sh                 # every dataset with stage-2 output
#   bash .../run_BENCHMARK_figures.sh pm_gt_persample # one dataset
# Datasets come from scripts/BENCHMARK_datasets.tsv. A dataset with no stage-2
# output is SKIPPED rather than failing, so the loop works while other tiers are
# still running.
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DELIV_DIR="$(dirname "$SCRIPT_DIR")"

# Activate the project env from conda's profile rather than ~/.bashrc: the
# interactive rc runs `module load gcc` and /etc/bashrc, both of which touch
# unbound variables and abort under `set -u`.
set +u
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
[ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ] && . "${CONDA_ROOT}/etc/profile.d/conda.sh"
[ -f "${CONDA_ROOT}/etc/profile.d/mamba.sh" ] && . "${CONDA_ROOT}/etc/profile.d/mamba.sh"
conda activate scRNAEditing
set -u

command -v Rscript >/dev/null || { echo "ERROR: Rscript not on PATH" >&2; exit 1; }

REG="${SCRIPT_DIR}/BENCHMARK_datasets.tsv"
[ -s "$REG" ] || { echo "ERROR: missing registry $REG" >&2; exit 1; }

if [ $# -ge 1 ]; then
    DATASETS="$*"
else
    DATASETS=$(awk -F'\t' 'NR>1 && $1 != "" {print $1}' "$REG")
fi

built=0
for DS in $DATASETS; do
    if [ ! -s "${DELIV_DIR}/data/${DS}/BENCHMARK_final_by_run.tsv.gz" ]; then
        echo "--- ${DS}: no stage-2 output, skipping"
        echo "      (build it with: Rscript ${SCRIPT_DIR}/BENCHMARK_prep_source_data.R ${DS})"
        continue
    fi
    echo "=== BENCHMARK_grid_source_data.R + BENCHMARK_grid_plot.R ${DS} ==="
    Rscript "${SCRIPT_DIR}/BENCHMARK_grid_source_data.R" "$DS"
    Rscript "${SCRIPT_DIR}/BENCHMARK_grid_plot.R" "$DS"
    built=$((built + 1))
done
[ "$built" -gt 0 ] || { echo "ERROR: no dataset had stage-2 output" >&2; exit 1; }

echo
echo "=== DONE ==="
echo "datasets built: ${built}"
echo "figures/      : $(find "${DELIV_DIR}/figures" -name '*.svg' 2>/dev/null | wc -l) SVG"
echo "source_data/  : $(find "${DELIV_DIR}/source_data" -name 'BENCHMARK_grid.*.tsv' ! -name '*.schema.tsv' 2>/dev/null | wc -l) TSV"

# figures/ is SVG-only by design; catch anything that slips back in.
stray=$(find "${DELIV_DIR}/figures" -type f ! -name '*.svg' | wc -l)
if [ "$stray" -gt 0 ]; then
    echo "WARNING: ${stray} non-SVG file(s) in figures/ — tables belong in source_data/" >&2
    find "${DELIV_DIR}/figures" -type f ! -name '*.svg' >&2
fi
