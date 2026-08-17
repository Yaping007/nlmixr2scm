#!/bin/bash
# ============================================================================
# submit_poweronly_est703_08132026.sh
# ----------------------------------------------------------------------------
# Orchestrator for the POWER-ONLY covariate-space rerun of nlmixr2-SCM.
#
# GOAL
#   Restrict the continuous-covariate shape menu to POWER ONLY (drop the
#   competing exponential/"lin" shape) and rerun the SCM benchmark:
#     N        = 80
#     scenarios= all 16
#     datasets = 1..100
#     structures = linCmt AND ode   (see WHY BOTH below)
#     estimator x optimizer = focei + bobyqa   (single active cell)
#     profile-on-stall warm-start = ON  (matches the production sweep)
#     shapes   = power   (--shapes "power")
#
# COMPETING vs POWER-ONLY -> SEPARATE OUTPUT TREE
#   The (est,opt) leaf dir is identical to the competing-shape sweep, so the two
#   WOULD collide in one tree. Per submit_one_array.sh guidance (pair a
#   non-default search knob with FORCE_RERUN=1 + a distinct OUT_ROOT), the
#   power-only results write to their OWN root and are aggregated separately:
#     output/scm_poweronly_est703_08132026
#
# WHY BOTH STRUCTURES
#   focei + bobyqa is a derivative-free outer optimiser, so (unlike lbfgsb3c)
#   linCmt does NOT degrade it. Both linCmt and ode are therefore informative
#   and are included. To run ONE structure only, edit STRUCTURES below
#   (e.g. STRUCTURES="ode").
#
# RESOURCES / THROTTLE
#   Default per-task LSF resources from bench_array.lsf (-n 4, -M 3000, -W 03:00)
#   and a %50 array throttle (100 datasets, 50 concurrent).
#
# USAGE (on HPCE, from repo root):
#   sed -i 's/\r$//' script/hpce_scm_estimator/*.sh script/hpce_scm_estimator/*.lsf
#   bash script/hpce_scm_estimator/submit_poweronly_est703_08132026.sh
#
#   # smoke-test one cell per structure first (recommended):
#   SMOKE=1 bash script/hpce_scm_estimator/submit_poweronly_est703_08132026.sh
# ============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- fixed configuration for this rerun ------------------------------------
# Edit STRUCTURES to "ode" or "linCmt" to run a single structure.
export STRUCTURES="linCmt ode"
export NS="80"
export SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
export FORCE_RERUN="1"                 # new (power-only) search -> regenerate all fits
export SHAPES="power"                  # <-- POWER-ONLY covariate space
export PROFILE_INIT_ON_STALL="TRUE"    # warm-start ON
NDS=100
MAXPAR=50

OUT_ROOT_POWERONLY="output/scm_poweronly_est703_08132026"

# ---- optional single-cell smoke test ---------------------------------------
if [ "${SMOKE:-0}" = "1" ]; then
  echo ">> SMOKE: one cell per structure (scn16, ds1) before the full sweep"
  for st in $STRUCTURES; do
    echo ">> SMOKE: focei + bobyqa , struct=$st , shapes=power -> $OUT_ROOT_POWERONLY"
    Rscript script/PerformanceEvaluation_scm_bench.R \
      --N 80 --scenario 16 --dataset 1 \
      --estimator focei --outer_opt bobyqa --structure "$st" \
      --out_root "$OUT_ROOT_POWERONLY" \
      --shapes "power" \
      --profile_init_on_stall TRUE --force_rerun
  done
  echo ">> SMOKE done. Inspect res_ds001.rds (selected rows should be power/cat only),"
  echo "   then re-run without SMOKE=1."
  exit 0
fi

# ---- full sweep: focei + bobyqa, warm ON, power-only -----------------------
echo ">> ===== POWER-ONLY (shapes=power, warm ON) -> $OUT_ROOT_POWERONLY ====="
OUT_ROOT="$OUT_ROOT_POWERONLY" \
  ESTIMATORS="focei" FOCEI_OPTS="bobyqa" \
  bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"

n_struct=$(echo $STRUCTURES | wc -w)
echo "=========================================================="
echo "Power-only SCM arms submitted."
echo "  root      = $OUT_ROOT_POWERONLY  (focei+bobyqa, warm ON, shapes=power)"
echo "  arrays    = 1 cell x $n_struct structure(s) x 16 scenarios"
echo "  each array = [1-$NDS]%$MAXPAR"
echo "Track:  bjobs -J 'bench_N80_scn*_*_focei_bobyqa'"
echo ""
echo "After completion, AGGREGATE this root:"
echo "  Rscript script/aggregate_scm_estimator2.1.R --root output \\"
echo "    --sub scm_poweronly_est703_08132026 \\"
echo "    --out_dir output/scm_poweronly_est703_08132026_aggregated"
echo "=========================================================="
