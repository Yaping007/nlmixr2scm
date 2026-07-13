#!/bin/bash
# ============================================================================
# submit_one_array.sh  (refit_estimators)
# ----------------------------------------------------------------------------
# Submit ONE LSF job array for one (N, scenario, estimator, outer_opt,
# boundary) cell of the TRUE-MODEL refit benchmark, array indices
# DS_START..DS_START+NDS-1.
#
# Usage:
#   bash script/hpce_refit_estimators/submit_one_array.sh \
#        <N> <scn> <est> <opt> <bnd> <n_ds> [maxpar] [ds_start]
#
# Example (scn16 N=80 foceif+lbfgsb3c tight, ds 1..10, 20 parallel):
#   bash script/hpce_refit_estimators/submit_one_array.sh 80 16 foceif lbfgsb3c tight 10 20
# ============================================================================

set -euo pipefail

SAMPLE_N=${1:?"Usage: $0 <N:40|80|300> <scn:1..16> <est> <opt|NA> <bnd> <n_ds> [maxpar] [ds_start]"}
SCN=${2:?"scenario required"}
EST=${3:?"estimator required (focei|foceif|irlsfoceif|saem)"}
OPT=${4:?"outer_opt required (bobyqa|nlminb|lbfgsb3c|NA)"}
BND=${5:?"boundary required (none|wide|narrow|tight)"}
NDS=${6:?"n_datasets required"}
MAXPAR=${7:-30}
DS_START=${8:-1}
DS_END=$((DS_START + NDS - 1))

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$HERE/../.." && pwd)}"

# Auto-detect scripts subdir holding the driver
SCRIPTS_DIR=""
for c in script scripts R; do
  if [ -f "$REPO_ROOT/$c/bench_refit_estimators.R" ]; then
    SCRIPTS_DIR="$c"; break
  fi
done
if [ -z "$SCRIPTS_DIR" ]; then
  echo "ERROR: bench_refit_estimators.R not found under $REPO_ROOT" >&2
  exit 1
fi

if [ ! -f "$HERE/refit_estimators_array.lsf" ]; then
  echo "ERROR: $HERE/refit_estimators_array.lsf missing" >&2
  exit 1
fi

SCN_PAD=$(printf '%02d' "$SCN")
JOBNAME="refitE_N${SAMPLE_N}_scn${SCN_PAD}_${EST}_${OPT}_${BND}"

LOGDIR="$REPO_ROOT/logs/refitE_N${SAMPLE_N}"
mkdir -p "$LOGDIR"

echo "Submitting: ${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}"
echo "  REPO_ROOT   = $REPO_ROOT"
echo "  SCRIPTS_DIR = $SCRIPTS_DIR"
echo "  logs        = $LOGDIR/${JOBNAME}.<jobid>.<idx>.{out,err}"

bsub \
  -J "${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}" \
  -o "${LOGDIR}/${JOBNAME}.%J.%I.out" \
  -e "${LOGDIR}/${JOBNAME}.%J.%I.err" \
  -env "all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, SAMPLE_N=${SAMPLE_N}, SCN=${SCN}, EST=${EST}, OPT=${OPT}, BND=${BND}" \
  < "$HERE/refit_estimators_array.lsf"
