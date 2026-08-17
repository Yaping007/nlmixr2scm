#!/bin/bash
# ============================================================================
# submit_lbfgsb3c_lincmt_est703_08152026.sh
# ----------------------------------------------------------------------------
# GROUND-TRUTH linCmt rerun of the gradient optimizer lbfgsb3c.
#
# WHY THIS EXISTS
#   The 5-arm ODE comparison already has all arms. The 5-arm *linCmt* (analytic
#   / Almquist) comparison is missing ONLY the two lbfgsb3c cells:
#       - focei + lbfgsb3c , linCmt , warm-up ON
#       - focei + lbfgsb3c , linCmt , warm-up OFF
#   bobyqa ON/OFF (linCmt) and PsN-SCM (advan4 = analytic) already exist.
#
#   NOTE ON COMPATIBILITY: lbfgsb3c is a gradient optimizer whose analytic
#   (Almquist) outer gradient only engages on an EXPLICIT ODE system. On the
#   analytic linCmt() solution the outer gradient falls back to finite
#   differences, so this arm is expected to be less well matched than on ODE.
#   We run it anyway to obtain the ground-truth linCmt behaviour.
#
# SEPARATE OUTPUT TREES (non-destructive)
#   These write to their OWN roots so they do NOT touch the existing ODE
#   lbfgsb3c roots/aggregates (which the ODE 5-arm figure depends on):
#     warm ON  -> output/lbfgsb3c_lincmt_est703_08152026_warmon
#     warm OFF -> output/lbfgsb3c_lincmt_est703_08152026_warmoff
#
# CONFIG (matches the ODE lbfgsb3c run)
#   N = 80, structure = linCmt only, all 16 scenarios x 100 datasets,
#   competing shapes (power + lin + cat) -- i.e. the DEFAULT search space
#   (do NOT set SHAPES; leave the competing menu on to match bobyqa/PsN arms).
#
# RESOURCES / THROTTLE
#   Default per-task LSF resources from bench_array.lsf (-n 4, -M 3000, -W 03:00)
#   and a %50 array throttle (100 datasets, 50 concurrent).
#
# USAGE (on HPCE, from repo root):
#   sed -i 's/\r$//' script/hpce_scm_estimator/*.sh script/hpce_scm_estimator/*.lsf
#   bash script/hpce_scm_estimator/submit_lbfgsb3c_lincmt_est703_08152026.sh
#
#   # smoke-test one cell per arm first (recommended):
#   SMOKE=1 bash script/hpce_scm_estimator/submit_lbfgsb3c_lincmt_est703_08152026.sh
# ============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- fixed configuration for this rerun ------------------------------------
export STRUCTURES="linCmt"             # <-- analytic structure ONLY
export NS="80"
export SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
export FORCE_RERUN="1"                 # new (linCmt) cells -> generate all fits
# NB: do NOT export SHAPES -> keep the default competing menu (power+lin+cat)
NDS=100
MAXPAR=50

ROOT_WARMON="output/lbfgsb3c_lincmt_est703_08152026_warmon"
ROOT_WARMOFF="output/lbfgsb3c_lincmt_est703_08152026_warmoff"

# ---- optional single-cell smoke test ---------------------------------------
if [ "${SMOKE:-0}" = "1" ]; then
  echo ">> SMOKE: one cell per arm (scn16, ds1, linCmt) before the full sweep"
  for warm in TRUE FALSE; do
    root=$([ "$warm" = "TRUE" ] && echo "$ROOT_WARMON" || echo "$ROOT_WARMOFF")
    echo ">> SMOKE: focei + lbfgsb3c , linCmt , warm=$warm -> $root"
    Rscript script/PerformanceEvaluation_scm_bench.R \
      --N 80 --scenario 16 --dataset 1 \
      --estimator focei --outer_opt lbfgsb3c --structure linCmt \
      --out_root "$root" \
      --profile_init_on_stall "$warm" --force_rerun
  done
  echo ">> SMOKE done. Inspect res_ds001.rds under both roots, then re-run without SMOKE=1."
  exit 0
fi

# ---- warm-start ON arm -----------------------------------------------------
echo ">> ===== linCmt WARM-START ON (profileInitOnStall = TRUE) -> $ROOT_WARMON ====="
OUT_ROOT="$ROOT_WARMON" PROFILE_INIT_ON_STALL="TRUE" \
  ESTIMATORS="focei" FOCEI_OPTS="lbfgsb3c" \
  bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"

# ---- warm-start OFF arm ----------------------------------------------------
echo ">> ===== linCmt WARM-START OFF (profileInitOnStall = FALSE) -> $ROOT_WARMOFF ====="
OUT_ROOT="$ROOT_WARMOFF" PROFILE_INIT_ON_STALL="FALSE" \
  ESTIMATORS="focei" FOCEI_OPTS="lbfgsb3c" \
  bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"

echo "=========================================================="
echo "linCmt lbfgsb3c arms submitted."
echo "  warm ON  -> $ROOT_WARMON   (focei+lbfgsb3c, linCmt)"
echo "  warm OFF -> $ROOT_WARMOFF  (focei+lbfgsb3c, linCmt)"
echo "  arrays   = 2 arms x 1 structure x 16 scenarios = 32 arrays"
echo "  each array = [1-$NDS]%$MAXPAR"
echo "Track:  bjobs -J 'bench_N80_scn*_linCmt_*lbfgsb3c*'"
echo ""
echo "After completion, AGGREGATE each root separately:"
echo "  Rscript script/aggregate_scm_estimator2.1.R --root output \\"
echo "    --sub lbfgsb3c_lincmt_est703_08152026_warmon \\"
echo "    --out_dir output/lbfgsb3c_lincmt_est703_08152026_warmon_aggregated"
echo "  Rscript script/aggregate_scm_estimator2.1.R --root output \\"
echo "    --sub lbfgsb3c_lincmt_est703_08152026_warmoff \\"
echo "    --out_dir output/lbfgsb3c_lincmt_est703_08152026_warmoff_aggregated"
echo "=========================================================="
