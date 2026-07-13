#!/bin/bash
# ============================================================================
# submit_all_arrays.sh
# ----------------------------------------------------------------------------
# Submit the estimator x optimizer bench grid.  Enumerates
# NS x SCENARIOS x valid (EST, OPT) cells; skips invalid combos.
#
# focei_bobyqa is preserved from the legacy migration (do NOT re-run); the
# focei row here defaults to {nlminb, lbfgsb3c} only. Override via
# FOCEI_OPTS='bobyqa nlminb lbfgsb3c' if you want to re-fit bobyqa too.
#
# VAE is not part of the grid (runSCM incompatibility on vae-returned ui).
#
# Env-var filters (space-separated lists):
#   NS='80'                          (default: 40 80 300)
#   SCENARIOS='16'                   (default: 1..16)
#   ESTIMATORS='focei foceif ...'    (default: focei foceif irlsfoceif saem)
#   FOCEI_OPTS='nlminb lbfgsb3c'     (default: nlminb lbfgsb3c)
#
# Usage:
#   bash script/hpce_bench/submit_all_arrays.sh <n_datasets> [maxpar]
#
# Example: pilot on scn16, N=80, all 6 new cells, 10 datasets:
#   NS=80 SCENARIOS=16 bash script/hpce_bench/submit_all_arrays.sh 10 20
#
# The 6 new cells submitted are:
#   focei_nlminb, focei_lbfgsb3c,
#   foceif_nlminb, foceif_lbfgsb3c,
#   irlsfoceif_lbfgsb3c, saem_NA
# (focei_bobyqa is preserved from migration; not resubmitted.)
# ============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_datasets> [maxpar]"}
MAXPAR=${2:-30}
DS_START=${DS_START:-1}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[ -n "${NS:-}"         ] || NS="40 80 300"
[ -n "${SCENARIOS:-}"  ] || SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
[ -n "${ESTIMATORS:-}" ] || ESTIMATORS="focei foceif saem"
[ -n "${FOCEI_OPTS:-}" ] || FOCEI_OPTS="nlminb lbfgsb3c"

# Valid outer_opts per estimator -- MUST mirror valid_combos() in
# estimator_factory.R.  focei row defaults to {nlminb, lbfgsb3c} to avoid
# re-running migrated focei_bobyqa; set FOCEI_OPTS to override.
# 2026-07-12: irlsfoceif excluded from default grid (no working cov step);
# opt into it by setting ESTIMATORS='irlsfoceif' explicitly.
declare -A OPT_FOR
OPT_FOR[focei]="$FOCEI_OPTS"
OPT_FOR[foceif]="nlminb lbfgsb3c"
OPT_FOR[irlsfoceif]="lbfgsb3c"
OPT_FOR[saem]="NA"

n_submitted=0
for N in $NS; do
  for SCN in $SCENARIOS; do
    for EST in $ESTIMATORS; do
      OPTS="${OPT_FOR[$EST]:-}"
      if [ -z "$OPTS" ]; then
        echo "WARN: no valid outer_opt for est=$EST -- skipping" >&2
        continue
      fi
      for OPT in $OPTS; do
        bash "$HERE/submit_one_array.sh" "$N" "$SCN" "$EST" "$OPT" "$NDS" "$MAXPAR" "$DS_START"
        n_submitted=$((n_submitted + 1))
      done
    done
  done
done

echo "=========================================================="
echo "Submitted $n_submitted arrays x $NDS datasets = $((n_submitted * NDS)) fits"
echo "Track:  bjobs -J 'bench_*'"
echo "Kill:   bkill  -J 'bench_*'"
echo "=========================================================="
