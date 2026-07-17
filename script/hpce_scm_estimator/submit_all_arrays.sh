#!/bin/bash
# ============================================================================
# submit_all_arrays.sh
# ----------------------------------------------------------------------------
# Submit the estimator x optimizer x structure bench grid.  Enumerates
# STRUCTURES x NS x SCENARIOS x valid (EST, OPT) cells.
# All cells now emit schema-2.1 records via package_scm_schema21() to
#   output/scm_bench/N{NN}/scn{SS}_{structure}/{est}_{opt}/res_ds{DDD}.rds
#
# Structure notes:
#   * linCmt : analytic 2-cmt oral. foceif/irlsfoceif DEGRADE to focei here
#              (no ODE to interaction-linearise); they are still run to
#              demonstrate 'prefer linCmt when possible'.
#   * ode    : explicit ODE; unlocks analytic gradients for foceif/irlsfoceif.
#
# focei_bobyqa is now part of the general sweep (FOCEI_OPTS includes bobyqa)
# so legacy linCmt focei_bobyqa is re-run under schema 2.1. The old
# output/full_scm_focei_bobyqa/ tree is left untouched.
#
# VAE is not part of this grid (runSCM incompatibility on vae-returned ui);
# VAE covariate selection has its own driver (vae_covsel_driver.R).
#
# Env-var filters (space-separated lists):
#   STRUCTURES='linCmt ode'          (default: linCmt ode)
#   NS='80'                          (default: 40 80 300)
#   SCENARIOS='16'                   (default: 1..16)
#   ESTIMATORS='focei foceif ...'    (default: focei foceif irlsfoceif)
#   FOCEI_OPTS='bobyqa nlminb ...'   (default: bobyqa nlminb lbfgsb3c)
#
# Usage:
#   bash script/hpce_scm_estimator/submit_all_arrays.sh <n_datasets> [maxpar]
#
# Example: pilot on scn16, N=80, ode only, 10 datasets:
#   STRUCTURES=ode NS=80 SCENARIOS=16 \
#     bash script/hpce_scm_estimator/submit_all_arrays.sh 10 20
# ============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_datasets> [maxpar]"}
MAXPAR=${2:-30}
DS_START=${DS_START:-1}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[ -n "${STRUCTURES:-}" ] || STRUCTURES="linCmt ode"
[ -n "${NS:-}"         ] || NS="40 80 300"
[ -n "${SCENARIOS:-}"  ] || SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
[ -n "${ESTIMATORS:-}" ] || ESTIMATORS="focei foceif irlsfoceif"
[ -n "${FOCEI_OPTS:-}" ] || FOCEI_OPTS="bobyqa nlminb lbfgsb3c"

# Valid outer_opts per estimator -- MUST mirror valid_combos() in
# estimator_factory.R (saem removed 2026-07-17).
declare -A OPT_FOR
OPT_FOR[focei]="$FOCEI_OPTS"
OPT_FOR[foceif]="nlminb lbfgsb3c"
OPT_FOR[irlsfoceif]="lbfgsb3c"

n_submitted=0
for STRUCTURE in $STRUCTURES; do
  for N in $NS; do
    for SCN in $SCENARIOS; do
      for EST in $ESTIMATORS; do
        OPTS="${OPT_FOR[$EST]:-}"
        if [ -z "$OPTS" ]; then
          echo "WARN: no valid outer_opt for est=$EST -- skipping" >&2
          continue
        fi
        for OPT in $OPTS; do
          bash "$HERE/submit_one_array.sh" "$N" "$SCN" "$EST" "$OPT" \
               "$STRUCTURE" "$NDS" "$MAXPAR" "$DS_START"
          n_submitted=$((n_submitted + 1))
        done
      done
    done
  done
done

echo "=========================================================="
echo "Submitted $n_submitted arrays x $NDS datasets = $((n_submitted * NDS)) fits"
echo "Track:  bjobs -J 'bench_*'"
echo "Kill:   bkill  -J 'bench_*'"
echo "=========================================================="
