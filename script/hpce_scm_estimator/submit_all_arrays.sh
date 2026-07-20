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
# Isolation knobs (for clean, contention-free timing):
#   CHAIN=1      run arrays strictly one-at-a-time (LSF ended-dependency chain)
#   EXCLUSIVE=1  request whole nodes (-x) so fork workers never starve
#
# Usage:
#   bash script/hpce_scm_estimator/submit_all_arrays.sh <n_datasets> [maxpar]
#
# Example: pilot on scn16, N=80, ode only, 10 datasets:
#   STRUCTURES=ode NS=80 SCENARIOS=16 \
#     bash script/hpce_scm_estimator/submit_all_arrays.sh 10 20
#
# Example: contention-free timing A/B (serial arrays, exclusive nodes):
#   CHAIN=1 EXCLUSIVE=1 STRUCTURES=ode NS=80 SCENARIOS=16 ESTIMATORS=focei \
#     FOCEI_OPTS=bobyqa \
#     bash script/hpce_scm_estimator/submit_all_arrays.sh 20 1
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

# --- serialization / isolation knobs ---------------------------------------
#   CHAIN=1     : submit arrays with an LSF dependency chain so only ONE array
#                 runs at a time (each waits for the previous to END). This
#                 eliminates cross-array host contention, giving clean wall
#                 times -- at the cost of throughput. Combine with a small
#                 MAXPAR (or EXCLUSIVE) so tasks WITHIN an array don't contend
#                 either.
#   EXCLUSIVE=1 : pass -x to each array so every task owns its whole node; the
#                 3 SCM fork workers are never starved -> trustworthy timing.
CHAIN=${CHAIN:-0}
export EXCLUSIVE=${EXCLUSIVE:-0}

if [ "$CHAIN" = "1" ]; then
  echo "CHAIN mode ON: arrays will run strictly sequentially (ended-dependency)."
fi
if [ "$EXCLUSIVE" = "1" ]; then
  echo "EXCLUSIVE mode ON: each array requests whole nodes (-x)."
fi

# temp file used to read back each array's job id for the dependency chain
JOBID_TMP="$(mktemp)"
trap 'rm -f "$JOBID_TMP"' EXIT
PREV_JOBID=""

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
          : > "$JOBID_TMP"
          if [ "$CHAIN" = "1" ] && [ -n "$PREV_JOBID" ]; then
            DEPEND_JOBID="$PREV_JOBID" JOBID_FILE="$JOBID_TMP" \
              bash "$HERE/submit_one_array.sh" "$N" "$SCN" "$EST" "$OPT" \
                   "$STRUCTURE" "$NDS" "$MAXPAR" "$DS_START"
          else
            JOBID_FILE="$JOBID_TMP" \
              bash "$HERE/submit_one_array.sh" "$N" "$SCN" "$EST" "$OPT" \
                   "$STRUCTURE" "$NDS" "$MAXPAR" "$DS_START"
          fi
          if [ "$CHAIN" = "1" ]; then
            PREV_JOBID="$(cat "$JOBID_TMP" 2>/dev/null || true)"
            [ -n "$PREV_JOBID" ] || \
              echo "WARN: could not read job id to chain on; next array will NOT wait." >&2
          fi
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
