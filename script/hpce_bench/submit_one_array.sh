#!/bin/bash
# ============================================================================
# submit_one_array.sh
# ----------------------------------------------------------------------------
# Submit ONE LSF job array for one (N, scenario, estimator, outer_opt) cell,
# with array indices 1..n_datasets.
#
# Usage:
#   bash script/hpce_bench/submit_one_array.sh <N> <scn> <est> <opt> <n_ds> [maxpar]
#
# Example (pilot: scn16, N=80, foceif+lbfgsb3c, 10 datasets, 20 parallel):
#   bash script/hpce_bench/submit_one_array.sh 80 16 foceif lbfgsb3c 10 20
# ============================================================================

set -euo pipefail

SAMPLE_N=${1:?"Usage: $0 <N:40|80|300> <scn:1..16> <est> <opt|NA> <n_ds> [maxpar] [ds_start]"}
SCN=${2:?"scenario required"}
EST=${3:?"estimator required (focei|foceif|irlsfoceif|saem|vae)"}
OPT=${4:?"outer_opt required (bobyqa|nlminb|lbfgsb3c|NA)"}
NDS=${5:?"n_datasets required"}
MAXPAR=${6:-30}
DS_START=${7:-1}
DS_END=$((DS_START + NDS - 1))

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$HERE/../.." && pwd)}"

# Auto-detect scripts subdir holding the driver
SCRIPTS_DIR=""
for c in script scripts R; do
  if [ -f "$REPO_ROOT/$c/PerformanceEvaluation_scm_bench.R" ]; then
    SCRIPTS_DIR="$c"; break
  fi
done
if [ -z "$SCRIPTS_DIR" ]; then
  echo "ERROR: PerformanceEvaluation_scm_bench.R not found under $REPO_ROOT" >&2
  exit 1
fi

if [ ! -f "$HERE/bench_array.lsf" ]; then
  echo "ERROR: $HERE/bench_array.lsf missing" >&2
  exit 1
fi

SCN_PAD=$(printf '%02d' "$SCN")
JOBNAME="bench_N${SAMPLE_N}_scn${SCN_PAD}_${EST}_${OPT}"

LOGDIR="$REPO_ROOT/logs/bench_N${SAMPLE_N}"
mkdir -p "$LOGDIR"

echo "Submitting: ${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}"
echo "  REPO_ROOT   = $REPO_ROOT"
echo "  SCRIPTS_DIR = $SCRIPTS_DIR"
echo "  logs        = $LOGDIR/${JOBNAME}.<jobid>.<idx>.{out,err}"

bsub \
  -J "${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}" \
  -o "${LOGDIR}/${JOBNAME}.%J.%I.out" \
  -e "${LOGDIR}/${JOBNAME}.%J.%I.err" \
  -env "all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, SAMPLE_N=${SAMPLE_N}, SCN=${SCN}, EST=${EST}, OPT=${OPT}" \
  < "$HERE/bench_array.lsf"
