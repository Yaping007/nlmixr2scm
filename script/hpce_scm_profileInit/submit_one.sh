#!/bin/bash
# ============================================================================
# submit_one.sh
# ----------------------------------------------------------------------------
# Submit ONE profileInit A/B task: a single (cell, arm, structure) on scn16.
#
# Usage:
#   bash script/hpce_profileInit/submit_one.sh <cell> <arm> <structure> [workers] [rxthreads] [cores] [walltime]
#
#   cell      : irlsfoceif_lbfgsb3c | focei_nlminb | focei_lbfgsb3c
#   arm       : profileOff | profileOn
#   structure : linCmt | ode
#   workers   : future multisession workers          (default 3)
#   rxthreads : rxode2 ODE threads per worker         (default 2)
#   cores     : #BSUB -n   (should be >= workers*rxthreads; default workers*rxthreads)
#   walltime  : #BSUB -W HH:MM                        (default 06:00)
#
# Examples:
#   # the heavy case you want to offload -- irlsfoceif on ODE, both arms:
#   bash script/hpce_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOff ode
#   bash script/hpce_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOn  ode
#
#   # request more muscle (6 workers x 2 threads = 12 cores, 10h):
#   bash script/hpce_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOn ode 6 2 12 10:00
# ============================================================================

set -euo pipefail

CELL=${1:?"Usage: $0 <cell> <arm> <structure> [workers] [rxthreads] [cores] [walltime]"}
ARM=${2:?"arm required (profileOff|profileOn)"}
STRUCTURE=${3:?"structure required (linCmt|ode)"}
WORKERS=${4:-3}
RXTHREADS=${5:-2}
CORES=${6:-$((WORKERS * RXTHREADS))}
WALLTIME=${7:-06:00}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$HERE/../.." && pwd)}"

# Auto-detect scripts subdir holding the A/B driver
SCRIPTS_DIR=""
for c in script scripts R; do
  if [ -f "$REPO_ROOT/$c/_ab_profileInit_scn16.R" ]; then
    SCRIPTS_DIR="$c"; break
  fi
done
if [ -z "$SCRIPTS_DIR" ]; then
  echo "ERROR: _ab_profileInit_scn16.R not found under $REPO_ROOT" >&2
  exit 1
fi
if [ ! -f "$HERE/profileInit_array.lsf" ]; then
  echo "ERROR: $HERE/profileInit_array.lsf missing" >&2
  exit 1
fi

JOBNAME="ab_scn16_${STRUCTURE}_${CELL}_${ARM}"
LOGDIR="$REPO_ROOT/logs/profileInit"
mkdir -p "$LOGDIR"

echo "Submitting: ${JOBNAME}"
echo "  cell/arm/struct = ${CELL} / ${ARM} / ${STRUCTURE}"
echo "  parallelism     = ${WORKERS} workers x ${RXTHREADS} rxThreads"
echo "  resources       = -n ${CORES}  -W ${WALLTIME}"
echo "  REPO_ROOT       = $REPO_ROOT"
echo "  SCRIPTS_DIR     = $SCRIPTS_DIR"
echo "  logs            = $LOGDIR/${JOBNAME}.<jobid>.{out,err}"

bsub \
  -J "${JOBNAME}" \
  -n "${CORES}" \
  -W "${WALLTIME}" \
  -o "${LOGDIR}/${JOBNAME}.%J.out" \
  -e "${LOGDIR}/${JOBNAME}.%J.err" \
  -env "all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, AB_CELL=${CELL}, AB_ARM=${ARM}, AB_STRUCTURE=${STRUCTURE}, AB_WORKERS=${WORKERS}, AB_RXTHREADS=${RXTHREADS}" \
  < "$HERE/profileInit_array.lsf"
