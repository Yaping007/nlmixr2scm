#!/bin/bash
# ============================================================================
# submit_lbfgsb3c_est703_08132026.sh
# ----------------------------------------------------------------------------
# One-shot orchestrator for the est-7.0.3 gradient-optimizer rerun.
#
# GOAL
#   N = 80, structure = ODE only, all 16 scenarios x 100 datasets, under the
#   manager-updated nlmixr2est (7.0.3 on HPCE), for TWO gradient estimators
#     (1) foceif + lbfgsb3c   -- analytic (Almquist) outer gradient
#     (2) focei  + lbfgsb3c   -- finite-difference outer gradient
#   each run BOTH with the 1-D profile-on-stall warm-start ON and OFF -> a
#   2 x 2 = FOUR-arm design:
#     A) foceif + lbfgsb3c , warm ON   (profileInitOnStall = TRUE)
#     B) focei  + lbfgsb3c , warm ON
#     C) foceif + lbfgsb3c , warm OFF  (profileInitOnStall = FALSE)
#     D) focei  + lbfgsb3c , warm OFF
#
# WARM ON vs OFF -> SEPARATE OUTPUT TREES
#   The (est,opt) leaf dir is identical for the warm-on and warm-off arms, so
#   they WOULD collide in one tree. Per submit_one_array.sh guidance ("pair
#   FALSE runs with FORCE_RERUN=1 + a distinct OUT_ROOT"), warm-on and warm-off
#   write to DISTINCT roots and are aggregated separately:
#     warm ON  -> output/lbfgsb3c_est703_08132026_warmon
#     warm OFF -> output/lbfgsb3c_est703_08132026_warmoff
#
# WHY ODE ONLY
#   lbfgsb3c is a gradient optimizer and the analytic gradient only engages on
#   an explicit ODE system; linCmt would degrade foceif to FD and make the two
#   estimators near-identical. ODE is the informative structure.
#
# RESOURCES / THROTTLE
#   Default per-task LSF resources from bench_array.lsf (-n 4, -M 3000, -W 03:00)
#   and a %50 array throttle (100 datasets, 50 concurrent).
#
# USAGE (on HPCE, from repo root):
#   sed -i 's/\r$//' script/hpce_scm_estimator/*.sh script/hpce_scm_estimator/*.lsf
#   bash script/hpce_scm_estimator/submit_lbfgsb3c_est703_08132026.sh
#
#   # smoke-test one cell per arm first (recommended):
#   SMOKE=1 bash script/hpce_scm_estimator/submit_lbfgsb3c_est703_08132026.sh
# ============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- fixed configuration for this rerun ------------------------------------
export STRUCTURES="ode"
export NS="80"
export SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
export FORCE_RERUN="1"                 # updated nlmixr2est -> regenerate all fits
NDS=100
MAXPAR=50

ROOT_WARMON="output/lbfgsb3c_est703_08132026_warmon"
ROOT_WARMOFF="output/lbfgsb3c_est703_08132026_warmoff"

# ---- optional single-cell smoke test ---------------------------------------
if [ "${SMOKE:-0}" = "1" ]; then
  echo ">> SMOKE: one cell per arm (scn16, ds1, ODE) before the full sweep"
  for warm in TRUE FALSE; do
    root=$([ "$warm" = "TRUE" ] && echo "$ROOT_WARMON" || echo "$ROOT_WARMOFF")
    for cfg in "foceif lbfgsb3c" "focei lbfgsb3c"; do
      set -- $cfg
      echo ">> SMOKE: $1 + $2 , warm=$warm -> $root"
      Rscript script/PerformanceEvaluation_scm_bench.R \
        --N 80 --scenario 16 --dataset 1 \
        --estimator "$1" --outer_opt "$2" --structure ode \
        --out_root "$root" \
        --profile_init_on_stall "$warm" --force_rerun
    done
  done
  echo ">> SMOKE done. Inspect res_ds001.rds under both roots, then re-run without SMOKE=1."
  exit 0
fi

# ---- warm-start ON arms (A, B) ---------------------------------------------
echo ">> ===== WARM-START ON (profileInitOnStall = TRUE) -> $ROOT_WARMON ====="
OUT_ROOT="$ROOT_WARMON" PROFILE_INIT_ON_STALL="TRUE" \
  ESTIMATORS="foceif" FOCEIF_OPTS="lbfgsb3c" \
  bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"
OUT_ROOT="$ROOT_WARMON" PROFILE_INIT_ON_STALL="TRUE" \
  ESTIMATORS="focei" FOCEI_OPTS="lbfgsb3c" \
  bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"

# ---- warm-start OFF arms (C, D) --------------------------------------------
echo ">> ===== WARM-START OFF (profileInitOnStall = FALSE) -> $ROOT_WARMOFF ====="
OUT_ROOT="$ROOT_WARMOFF" PROFILE_INIT_ON_STALL="FALSE" \
  ESTIMATORS="foceif" FOCEIF_OPTS="lbfgsb3c" \
  bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"
OUT_ROOT="$ROOT_WARMOFF" PROFILE_INIT_ON_STALL="FALSE" \
  ESTIMATORS="focei" FOCEI_OPTS="lbfgsb3c" \
  bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"

echo "=========================================================="
echo "All four arms submitted."
echo "  warm ON  -> $ROOT_WARMON   (foceif+lbfgsb3c, focei+lbfgsb3c)"
echo "  warm OFF -> $ROOT_WARMOFF  (foceif+lbfgsb3c, focei+lbfgsb3c)"
echo "  arrays   = 4 arms x 1 structure x 16 scenarios = 64 arrays"
echo "  each array = [1-$NDS]%$MAXPAR"
echo "Track:  bjobs -J 'bench_N80_scn*_ode_*lbfgsb3c*'"
echo ""
echo "After completion, AGGREGATE each root separately:"
echo "  Rscript script/aggregate_scm_estimator2.1.R --root output \\"
echo "    --sub lbfgsb3c_est703_08132026_warmon \\"
echo "    --out_dir output/lbfgsb3c_est703_08132026_warmon_aggregated"
echo "  Rscript script/aggregate_scm_estimator2.1.R --root output \\"
echo "    --sub lbfgsb3c_est703_08132026_warmoff \\"
echo "    --out_dir output/lbfgsb3c_est703_08132026_warmoff_aggregated"
echo ""
echo "FD-switch trace (analytic<->FD fallback), one dataset per estimator:"
echo "  DIAG_STRUCTURE=ode Rscript script/diagnose_scm_gradient_trace.R"
echo "=========================================================="
