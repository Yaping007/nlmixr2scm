#!/bin/bash
# ============================================================================
# submit_timing_ab.sh
# ----------------------------------------------------------------------------
# Fork-parallelism timing A/B: run identical (N, scenario, structure) cells
# with workers=1 (serial SCM) vs workers=3 (forked SCM), on EXCLUSIVE nodes,
# so the ONLY difference is the fork width and the wall time is uncontended.
#
# The publishable metric is the SCM-phase wall-time ratio:
#     speedup = scm_sec(workers=1) / scm_sec(workers=3)
# read from output/scm_timing_ab/w1 vs .../w3 via aggregate_scm_timing.R.
#
# Each arm goes to its OWN output tree + job-name tag so tier-1 caching and
# LSF logs never collide:
#     arm w1 -> OUT_ROOT=output/scm_timing_ab/w1  JOBTAG=w1  WORKERS=1
#     arm w3 -> OUT_ROOT=output/scm_timing_ab/w3  JOBTAG=w3  WORKERS=3
# Both pin RX_THREADS=1 (fork-safe; keeps the comparison about workers only).
#
# Isolation: EXCLUSIVE=1 (whole node) + MAXPAR=1 (one task at a time) + the
# arms are chained (w3 waits for w1) so nothing contends.
#
# Usage:
#   bash script/hpce_scm_estimator/submit_timing_ab.sh <n_datasets>
#
# Env filters (space-separated; small grids recommended for the A/B):
#   STRUCTURES (default 'linCmt ode')  NS (default '40 80 300')
#   SCENARIOS  (default 16)            EST (default focei)  OPT (default bobyqa)
#
# Example (default focei_bobyqa, scn16, all N & structures, 20 datasets):
#   bash script/hpce_scm_estimator/submit_timing_ab.sh 20
# ============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_datasets>"}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

STRUCTURES=${STRUCTURES:-"linCmt ode"}
NS=${NS:-"40 80 300"}
SCENARIOS=${SCENARIOS:-16}
EST=${EST:-focei}
OPT=${OPT:-bobyqa}
# Concurrency per arm. Safe with EXCLUSIVE=1 because each task owns its own
# node -> %MAXPAR tasks land on MAXPAR different nodes and cannot contend.
MAXPAR=${MAXPAR:-1}

echo "Timing A/B: EST=$EST OPT=$OPT  N={$NS}  scn={$SCENARIOS}  struct={$STRUCTURES}  nds=$NDS"
echo "Arms: w1 (serial) vs w3 (fork x3); exclusive nodes, serial, chained."

# Run each arm as its own serialized+exclusive sweep. CHAIN inside each arm
# serializes its cells; we also chain arm w3 after arm w1 by waiting on w1 to
# drain (simple approach: submit w1 fully, then w3 -- both are already serial
# and exclusive, so cross-arm overlap is limited and harmless for timing since
# each task still owns its node).
submit_arm() {
  local workers="$1" tag="$2"
  echo "==== arm ${tag}: WORKERS=${workers} ===="
  CHAIN=1 EXCLUSIVE=1 \
  WORKERS="$workers" RX_THREADS=1 JOBTAG="$tag" \
  OUT_ROOT="output/scm_timing_ab/${tag}" \
  STRUCTURES="$STRUCTURES" NS="$NS" SCENARIOS="$SCENARIOS" \
  ESTIMATORS="$EST" FOCEI_OPTS="$OPT" \
    bash "$HERE/submit_all_arrays.sh" "$NDS" "$MAXPAR"
}

submit_arm 1 w1
submit_arm 3 w3

echo "=========================================================="
echo "Submitted both arms. When done, extract the speedup:"
echo "  Rscript script/aggregate_scm_timing.R --root output/scm_timing_ab/w1 --out_dir output/scm_timing_ab/agg_w1"
echo "  Rscript script/aggregate_scm_timing.R --root output/scm_timing_ab/w3 --out_dir output/scm_timing_ab/agg_w3"
echo "  # then compare wall_scm_med (w1) / wall_scm_med (w3) per (N, structure)"
echo "=========================================================="
