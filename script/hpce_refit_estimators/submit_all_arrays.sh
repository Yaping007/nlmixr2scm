#!/bin/bash
# ============================================================================
# submit_all_arrays.sh  (refit_estimators)
# ----------------------------------------------------------------------------
# Submit the TRUE-MODEL refit grid across estimators.  Enumerates
# NS x SCENARIOS x BOUNDARIES x valid (EST, OPT); skips invalid combos.
#
# Env-var filters (space-separated lists):
#   NS='80'                             (default: 40 80 300)
#   SCENARIOS='16'                      (default: 1..16)
#   BOUNDARIES='tight'                  (default: tight)
#   ESTIMATORS='focei foceif ...'       (default: focei foceif irlsfoceif saem)
#   FOCEI_OPTS='nlminb lbfgsb3c bobyqa' (default: nlminb lbfgsb3c bobyqa)
#   DS_START=1                          (starting dataset index; default 1)
#
# Usage:
#   bash script/hpce_refit_estimators/submit_all_arrays.sh <n_datasets> [maxpar]
#
# Example (scn16 N=80 tight, ds 1..10, 20 parallel, all 7 cells):
#   NS=80 SCENARIOS=16 BOUNDARIES=tight \
#     bash script/hpce_refit_estimators/submit_all_arrays.sh 10 20
#
# Cells submitted (valid_combos() in estimator_factory.R):
#   focei_bobyqa, focei_nlminb, focei_lbfgsb3c,
#   foceif_nlminb, foceif_lbfgsb3c,
#   irlsfoceif_lbfgsb3c, saem_NA
# ============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_datasets> [maxpar]"}
MAXPAR=${2:-30}
DS_START=${DS_START:-1}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[ -n "${NS:-}"         ] || NS="40 80 300"
[ -n "${SCENARIOS:-}"  ] || SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
[ -n "${BOUNDARIES:-}" ] || BOUNDARIES="tight"
[ -n "${ESTIMATORS:-}" ] || ESTIMATORS="focei foceif saem"
[ -n "${FOCEI_OPTS:-}" ] || FOCEI_OPTS="nlminb lbfgsb3c bobyqa"

declare -A OPT_FOR
OPT_FOR[focei]="$FOCEI_OPTS"
OPT_FOR[foceif]="nlminb lbfgsb3c"
OPT_FOR[irlsfoceif]="lbfgsb3c"
OPT_FOR[saem]="NA"

n_submitted=0
for N in $NS; do
  for SCN in $SCENARIOS; do
    for BND in $BOUNDARIES; do
      for EST in $ESTIMATORS; do
        OPTS="${OPT_FOR[$EST]:-}"
        if [ -z "$OPTS" ]; then
          echo "WARN: no valid outer_opt for est=$EST -- skipping" >&2
          continue
        fi
        for OPT in $OPTS; do
          bash "$HERE/submit_one_array.sh" \
               "$N" "$SCN" "$EST" "$OPT" "$BND" "$NDS" "$MAXPAR" "$DS_START"
          n_submitted=$((n_submitted + 1))
        done
      done
    done
  done
done

echo "=========================================================="
echo "Submitted $n_submitted arrays x $NDS datasets = $((n_submitted * NDS)) fits"
echo "Track:  bjobs -J 'refitE_*'"
echo "Kill:   bkill  -J 'refitE_*'"
echo "=========================================================="
