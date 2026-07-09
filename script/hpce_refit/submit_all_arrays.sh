#!/bin/bash
# ==============================================================================
# submit_all_arrays.sh
# ------------------------------------------------------------------------------
# Submit ALL LSF job arrays for the refit-true-model study.  One array per
# (cohort, scenario, boundary) triple:
#
#   scenario 1        -> skip (no covariates -> no thetas to bound)
#   scenario 2        -> single boundary ("none"): only categorical TH_SEX_VC,
#                        which is always unbounded regardless of config
#   scenarios 3..16   -> three boundaries: none, wide, narrow
#
# Combos per cohort: 1 + 14 * 3 = 43
# Cohorts:           N40, N80, N300  (3)
# Total arrays:      43 * 3 = 129
# Total fits (pilot, 10 ds):  129 * 10  = 1,290
# Total fits (full,  250 ds): 129 * 250 = 32,250
#
# Usage:
#   bash <hpce_dir>/submit_all_arrays.sh <n_datasets> [max_parallel] [cohorts...]
#
# The script auto-detects its own directory, so <hpce_dir> can be
# scripts/hpce/, script/hpce_refit/, or any path you choose.
# Override REPO_ROOT via env var if the scripts are not exactly 2 levels
# below the repo root.
#
# Examples:
#   # Pilot (10 datasets per combo, all 3 cohorts):
#   bash script/hpce_refit/submit_all_arrays.sh 10 50
#
#   # Full run (250 datasets per combo, all 3 cohorts):
#   bash script/hpce_refit/submit_all_arrays.sh 250 100
#
#   # Only N=300 cohort, 250 datasets:
#   bash script/hpce_refit/submit_all_arrays.sh 250 100 N300
# ==============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_datasets> [max_parallel] [cohorts...]  (10 pilot, 250 full)"}
MAXPAR=${2:-50}
shift 2 2>/dev/null || shift $#

# Remaining args = cohort list; default to all three
if [ $# -gt 0 ]; then
    COHORTS=("$@")
else
    COHORTS=(N40 N80 N300)
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- (scenario, boundaries) map --------------------------------------------
# Use two parallel arrays (bash 3 compatible; assoc arrays need bash 4+)
SCENARIOS=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
BOUNDARIES_FOR_2="none"
BOUNDARIES_FOR_MULTI="none wide narrow"

# ---- Sanity check -----------------------------------------------------------
if ! [[ "$NDS" =~ ^[0-9]+$ ]] || [ "$NDS" -lt 1 ]; then
    echo "ERROR: n_datasets must be a positive integer, got '$NDS'" >&2
    exit 1
fi

echo "=========================================================="
echo "Submitting refit-true-model arrays"
echo "  n_datasets  : $NDS"
echo "  max_parallel: $MAXPAR"
echo "  cohorts     : ${COHORTS[*]}"
echo "  Total combos: $(( ${#COHORTS[@]} * 43 ))"
echo "  Total fits  : $(( ${#COHORTS[@]} * 43 * NDS ))"
echo "=========================================================="

n_submitted=0
for COHORT in "${COHORTS[@]}"; do
    for SCN in "${SCENARIOS[@]}"; do
        if [ "$SCN" -eq 2 ]; then
            BOUNDARIES="$BOUNDARIES_FOR_2"
        else
            BOUNDARIES="$BOUNDARIES_FOR_MULTI"
        fi
        for BND in $BOUNDARIES; do
            bash "$HERE/submit_one_array.sh" "$COHORT" "$SCN" "$BND" "$NDS" "$MAXPAR"
            n_submitted=$((n_submitted + 1))
        done
    done
done

echo "=========================================================="
echo "Submitted $n_submitted arrays."
echo "Track with:  bjobs -J 'refit_*'"
echo "Kill all with: bkill -J 'refit_*'"
echo "=========================================================="
