#!/bin/bash
# ==============================================================================
# submit_all_arrays.sh
# ------------------------------------------------------------------------------
# Submit ALL LSF job arrays for the refit-true-model study.  One array per
# (cohort, scenario, boundary) triple.
#
# Boundary tiers (see script/true_model_factory.R):
#   * none    -> no bounds on continuous cov thetas
#   * wide    -> c(-1e5, 1e5)       PsN default
#   * narrow  -> c(-10,  10)        paper's narrow tier
#   * tight   -> c(-5,   5)         sensitivity tier (added 2026-07-09)
#
# Scenario rules:
#   scenario 1      -> skipped entirely (no covariates -> no thetas to bound)
#   scenario 2      -> only categorical TH_SEX_VC (always unbounded); the only
#                      meaningful boundary is "none". Submitted only when the
#                      filter includes "none".
#   scenarios 3..16 -> honor whatever is in BOUNDARIES (default: none wide narrow)
#
# Filter via ENV VARS (preferred over positional cohort args):
#   COHORTS='N300'          -> only submit that cohort (default: N40 N80 N300)
#   BOUNDARIES='narrow'     -> only that boundary  (default: none wide narrow)
#   BOUNDARIES='none wide'  -> multiple; space-separated
#   BOUNDARIES='tight'      -> the new sensitivity tier only
#
# The script auto-detects its own directory, so <hpce_dir> can be
# scripts/hpce/, script/hpce_refit/, or any path you choose.
# Override REPO_ROOT via env var if the scripts are not exactly 2 levels
# below the repo root.
#
# Usage:
#   bash <hpce_dir>/submit_all_arrays.sh <n_datasets> [max_parallel] [cohorts...]
#
# Examples:
#   # Pilot (10 datasets per combo, all 3 cohorts, paper's 3 boundaries):
#   bash script/hpce_refit/submit_all_arrays.sh 10 50
#
#   # Full run (250 datasets per combo, all 3 cohorts, paper's 3 boundaries):
#   bash script/hpce_refit/submit_all_arrays.sh 250 100
#
#   # Only N=300 cohort, 250 datasets (positional; all 3 boundaries):
#   bash script/hpce_refit/submit_all_arrays.sh 250 100 N300
#
#   # Only N=300 x narrow, 250 datasets (env-var filter):
#   BOUNDARIES=narrow COHORTS=N300 bash script/hpce_refit/submit_all_arrays.sh 250 50
#
#   # Only N=300 x tight (new sensitivity tier), 250 datasets:
#   BOUNDARIES=tight  COHORTS=N300 bash script/hpce_refit/submit_all_arrays.sh 250 50
# ==============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_datasets> [max_parallel] [cohorts...]  (10 pilot, 250 full)"}
MAXPAR=${2:-50}
shift 2 2>/dev/null || shift $#

# ---- Cohort list resolution ------------------------------------------------
# Priority:
#   1. Positional args (backward-compat)
#   2. COHORTS env var (space-separated)
#   3. Default: all three
if [ $# -gt 0 ]; then
    COHORTS_ARR=("$@")
elif [ -n "${COHORTS:-}" ]; then
    read -r -a COHORTS_ARR <<< "$COHORTS"
else
    COHORTS_ARR=(N40 N80 N300)
fi

# ---- Boundary filter (env-var only) ---------------------------------------
# Default = paper's 3 tiers (none wide narrow).  The new "tight" tier
# is NOT in the default set -- request it explicitly via BOUNDARIES=tight.
if [ -n "${BOUNDARIES:-}" ]; then
    read -r -a BND_FILTER <<< "$BOUNDARIES"
else
    BND_FILTER=(none wide narrow)
fi

# Validate every requested boundary against the master list.
ALL_BOUNDARIES=(none wide narrow tight)
_is_valid_bnd() {
    local needle="$1" b
    for b in "${ALL_BOUNDARIES[@]}"; do
        [ "$b" = "$needle" ] && return 0
    done
    return 1
}
for b in "${BND_FILTER[@]}"; do
    if ! _is_valid_bnd "$b"; then
        echo "ERROR: boundary '$b' not one of: ${ALL_BOUNDARIES[*]}" >&2
        exit 1
    fi
done

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Scenarios --------------------------------------------------------------
SCENARIOS=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

# ---- Sanity check -----------------------------------------------------------
if ! [[ "$NDS" =~ ^[0-9]+$ ]] || [ "$NDS" -lt 1 ]; then
    echo "ERROR: n_datasets must be a positive integer, got '$NDS'" >&2
    exit 1
fi

# ---- Count expected submissions for the banner ----------------------------
# Rules:
#   scn == 2 -> only submit when "none" is in the filter
#   scn > 2  -> submit for every boundary in the filter
n_expected=0
for SCN in "${SCENARIOS[@]}"; do
    for BND in "${BND_FILTER[@]}"; do
        if [ "$SCN" -eq 2 ] && [ "$BND" != "none" ]; then
            continue
        fi
        n_expected=$((n_expected + ${#COHORTS_ARR[@]}))
    done
done

echo "=========================================================="
echo "Submitting refit-true-model arrays"
echo "  n_datasets  : $NDS"
echo "  max_parallel: $MAXPAR"
echo "  cohorts     : ${COHORTS_ARR[*]}"
echo "  boundaries  : ${BND_FILTER[*]}"
echo "  arrays      : $n_expected"
echo "  total fits  : $((n_expected * NDS))"
echo "=========================================================="

if [ "$n_expected" -eq 0 ]; then
    echo "WARNING: 0 arrays match the filter. Nothing to submit." >&2
    exit 0
fi

n_submitted=0
for COHORT in "${COHORTS_ARR[@]}"; do
    for SCN in "${SCENARIOS[@]}"; do
        for BND in "${BND_FILTER[@]}"; do
            if [ "$SCN" -eq 2 ] && [ "$BND" != "none" ]; then
                continue
            fi
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
