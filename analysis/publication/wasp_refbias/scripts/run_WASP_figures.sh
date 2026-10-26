#!/bin/bash
# ============================================================================
# run_WASP_figures.sh — rebuild every figure in the WASP reference-bias
# supplement from the prepared source data, for one dataset or all of them.
# ============================================================================
# Assumes the source-data chain has already been run for that dataset:
#   0.  bash/qsub human/<ds>/scripts/submit_wasp_realign_*.sh  (STAR-WASP)
#   0b. qsub    scripts/WASP_filter_bams.sh    <ds>   (vW filter -> passed BAMs)
#   1.  qsub    scripts/WASP_prep_pileups.sh   <ds>   (BAMs -> pileups)
#   2.  Rscript scripts/WASP_prep_source_data.R <ds>  (pileups -> tables)
# This script is step 3 and is cheap (seconds per panel), so it is safe to run
# interactively on a login node.
#
# WASP_refbias_source_data.R writes the plotted tables and WASP_refbias_plot.R
# draws each lettered panel and the combined figure they assemble into from
# those tables only. Every SVG goes to figures/<ds>/; every TSV (source data,
# .sources.tsv sidecars and figure_source_map.tsv) goes to source_data/<ds>/.
#
# Usage:
#   bash .../run_WASP_figures.sh            # every dataset with source data
#   bash .../run_WASP_figures.sh jzlss2     # just one
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

# With no argument, build every dataset whose stage-2 output exists. Datasets
# in the registry that have not been run yet are skipped, not failed: the
# registry lists what the diagnostic CAN cover, which is deliberately ahead of
# what has been aligned.
if [ $# -ge 1 ]; then
    DATASETS=("$@")
else
    mapfile -t DATASETS < <(
        awk -F'\t' 'NR>1 {print $1}' "${SCRIPT_DIR}/WASP_datasets.tsv" |
        while read -r ds; do
            [ -s "${DELIV_DIR}/data/${ds}/WASP_persite_evaluable.tsv.gz" ] && echo "$ds"
        done
    )
    [ ${#DATASETS[@]} -ge 1 ] || {
        echo "ERROR: no dataset has stage-2 output under ${DELIV_DIR}/data/" >&2
        echo "  Run: Rscript ${SCRIPT_DIR}/WASP_prep_source_data.R <dataset>" >&2
        exit 1
    }
fi

for ds in "${DATASETS[@]}"; do
    if [ ! -s "${DELIV_DIR}/data/${ds}/WASP_persite_evaluable.tsv.gz" ]; then
        echo "ERROR: source data missing for ${ds}. Run first:" >&2
        echo "  Rscript ${SCRIPT_DIR}/WASP_prep_source_data.R ${ds}" >&2
        exit 1
    fi
    echo "=== WASP_refbias_source_data.R + WASP_refbias_plot.R ${ds} ==="
    Rscript "${SCRIPT_DIR}/WASP_refbias_source_data.R" "$ds"
    Rscript "${SCRIPT_DIR}/WASP_refbias_plot.R" "$ds"
    echo
done

# Cross-platform summary figure. Needs >=2 datasets AND the site-clustered
# intervals from WASP_gap_ci.sh, which is an SGE job — so this is skipped with
# a notice rather than failing the whole rebuild.
if [ ${#DATASETS[@]} -ge 2 ]; then
    missing=""
    for ds in "${DATASETS[@]}"; do
        [ -s "${DELIV_DIR}/report/${ds}/WASP_gap_ci.tsv" ] || missing="${missing} ${ds}"
    done
    if [ -n "$missing" ]; then
        echo "SKIP cross-platform figure: no WASP_gap_ci.tsv for${missing}" >&2
        echo "  qsub ${SCRIPT_DIR}/WASP_gap_ci.sh${missing}" >&2
    else
        echo "=== WASP_refbias_platforms_source_data.R ${DATASETS[*]} ==="
        Rscript "${SCRIPT_DIR}/WASP_refbias_platforms_source_data.R" "${DATASETS[@]}"
        echo "=== WASP_refbias_platforms_plot.R ==="
        Rscript "${SCRIPT_DIR}/WASP_refbias_platforms_plot.R"
        echo
    fi
fi

echo "=== DONE — ${#DATASETS[@]} dataset(s): ${DATASETS[*]} ==="
for ds in "${DATASETS[@]}"; do
    printf '  %-8s %s SVG | %s source-data TSV\n' "$ds" \
        "$(ls -1 "${DELIV_DIR}/figures/${ds}"/*.svg 2>/dev/null | wc -l)" \
        "$(ls -1 "${DELIV_DIR}/source_data/${ds}"/WASP_refbias.*.tsv* 2>/dev/null | grep -vc '\.schema\.tsv$')"
done

# figures/ is SVG-only by design; catch anything that slips back in.
stray=$(find "${DELIV_DIR}/figures" -type f ! -name '*.svg' | wc -l)
if [ "$stray" -gt 0 ]; then
    echo "WARNING: ${stray} non-SVG file(s) under figures/ — tables belong in source_data/" >&2
    find "${DELIV_DIR}/figures" -type f ! -name '*.svg' >&2
fi
