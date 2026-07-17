#!/bin/bash
# ============================================================================
# submit_all_arrays.sh
# ----------------------------------------------------------------------------
# Enumerate the VAE covariate-selection pilot grid and submit one LSF array
# per (N, scenario, structure) cell.
#
#   Grid : NS x SCENARIOS x STRUCTURES  (one array each; index = dataset)
#   Full pilot default: 3 N x 16 scn x 2 struct = 96 arrays, x NDS datasets.
#
# Usage:
#   bash script/hpce_vae_covsel/submit_all_arrays.sh [n_ds] [maxpar] [ds_start]
#
#   n_ds     datasets per cell   (default 5   -- the pilot size)
#   maxpar   array throttle %N   (default 20)
#   ds_start first dataset index (default 1)
#
# Override any axis via env vars (space-separated):
#   NS="40 80 300"   SCENARIOS="1 .. 16"   STRUCTURES="linCmt ode"
#
# Examples
#   # SINGLE PROBE first (recommended): N80, scn16, linCmt, 5 datasets
#   NS=80 SCENARIOS=16 STRUCTURES=linCmt \
#     bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
#
#   # FULL pilot (only after the probe's resources look OK):
#   bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
# ============================================================================

set -euo pipefail

NDS=${1:-5}
MAXPAR=${2:-20}
DS_START=${3:-1}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[ -n "${NS:-}" ]         || NS="40 80 300"
[ -n "${SCENARIOS:-}" ]  || SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
[ -n "${STRUCTURES:-}" ] || STRUCTURES="linCmt ode"

echo "=========================================================="
echo "VAE covsel pilot submission"
echo "  NS         = $NS"
echo "  SCENARIOS  = $SCENARIOS"
echo "  STRUCTURES = $STRUCTURES"
echo "  NDS        = $NDS   (datasets per cell)"
echo "  MAXPAR     = $MAXPAR"
echo "  DS_START   = $DS_START"
echo "=========================================================="

n_cells=0
for N in $NS; do
  for SCN in $SCENARIOS; do
    for ST in $STRUCTURES; do
      bash "$HERE/submit_one_array.sh" "$N" "$SCN" "$ST" "$NDS" "$MAXPAR" "$DS_START"
      n_cells=$((n_cells + 1))
    done
  done
done

echo "=========================================================="
echo "Submitted $n_cells array(s), $((n_cells * NDS)) task(s) total."
echo "Monitor: bjobs -J 'vaecov_*'"
echo "=========================================================="
