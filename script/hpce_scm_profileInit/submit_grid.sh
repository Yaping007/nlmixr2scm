#!/bin/bash
# ============================================================================
# submit_grid.sh
# ----------------------------------------------------------------------------
# Fan out the whole profileInit A/B grid: every (cell x arm) for one structure,
# each as an independent LSF job (so they run concurrently -- the whole point
# of offloading to HPCE).
#
# Usage:
#   bash script/hpce_profileInit/submit_grid.sh <structure> [workers] [rxthreads] [cores] [walltime]
#
# Examples:
#   # the run that was too slow locally -- all 3 cells x 2 arms on ODE:
#   bash script/hpce_profileInit/submit_grid.sh ode
#
#   # linCmt grid with beefier tasks:
#   bash script/hpce_profileInit/submit_grid.sh linCmt 6 2 12 04:00
#
# Cells and arms are the same as the A/B driver's `cells` table.
# ============================================================================

set -euo pipefail

STRUCTURE=${1:?"Usage: $0 <structure:linCmt|ode> [workers] [rxthreads] [cores] [walltime]"}
WORKERS=${2:-3}
RXTHREADS=${3:-2}
CORES=${4:-$((WORKERS * RXTHREADS))}
WALLTIME=${5:-06:00}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CELLS=(irlsfoceif_lbfgsb3c focei_nlminb focei_lbfgsb3c)
ARMS=(profileOff profileOn)

echo "Submitting profileInit A/B grid: structure=${STRUCTURE}, "
echo "  ${#CELLS[@]} cells x ${#ARMS[@]} arms = $(( ${#CELLS[@]} * ${#ARMS[@]} )) jobs"
echo "  each: ${WORKERS}w x ${RXTHREADS}t (-n ${CORES}, -W ${WALLTIME})"
echo

for cell in "${CELLS[@]}"; do
  for arm in "${ARMS[@]}"; do
    bash "$HERE/submit_one.sh" "$cell" "$arm" "$STRUCTURE" \
         "$WORKERS" "$RXTHREADS" "$CORES" "$WALLTIME"
    echo
  done
done

echo "All jobs submitted. Track with:  bjobs -w"
