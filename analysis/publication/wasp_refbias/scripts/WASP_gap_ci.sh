#!/bin/bash
# ============================================================================
# WASP_gap_ci.sh — SGE wrapper for WASP_gap_ci.R
# ============================================================================
# >>> DO NOT RUN THE R SCRIPT ON THE LOGIN NODE. <<<
# It reads ~1 GB of pileup TSV per dataset and bootstraps 500 resamples over
# ~1.8 M sites, twice. That is minutes of saturated CPU and several GB of RAM
# per dataset — cluster work, not login-node work.
#
# Usage:  qsub publication/wasp_refbias/scripts/WASP_gap_ci.sh [dataset ...]
#         (no argument = jzlxv3 and jzlss2)
#
#$ -N WASP_gap_ci
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 4
#$ -l mem=48G
#$ -l h_rt=8:00:00
# `-l mem`, not `-l h_vmem`: this is a hand-written qsub, so the config sweep
# that fixed the scheduler clamp (see CHANGELOG) does not cover it.
# Do NOT pass -V (SGE JSV / truncated PATH).

set -euo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

# SGE copies the job script into a spool dir, so ${BASH_SOURCE[0]} does not
# point into the repo. Prefer the real location interactively, fall back under
# qsub. (Same trap that killed job 1710384.)
WASP_SCRIPTS_DEFAULT="${SCRNA_DIR:-.}/publication/wasp_refbias/scripts"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/WASP_gap_ci.R" ] || SCRIPT_DIR="${WASP_SCRIPTS_DEFAULT}"
[ -f "${SCRIPT_DIR}/WASP_gap_ci.R" ] || {
    echo "ERROR: cannot locate WASP_gap_ci.R (tried ${SCRIPT_DIR})" >&2; exit 1; }

echo "=== $(date) === WASP gap CI (site-clustered) ==="
Rscript "${SCRIPT_DIR}/WASP_gap_ci.R" "$@"
echo "=== $(date) === DONE ==="
