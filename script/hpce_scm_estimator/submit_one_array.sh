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

SAMPLE_N=${1:?"Usage: $0 <N:40|80|300> <scn:1..16> <est> <opt|NA> <structure:linCmt|ode> <n_ds> [maxpar] [ds_start]"}
SCN=${2:?"scenario required"}
EST=${3:?"estimator required (focei|foceif|irlsfoceif|vae)"}
OPT=${4:?"outer_opt required (bobyqa|nlminb|lbfgsb3c|NA)"}
STRUCTURE=${5:?"structure required (linCmt|ode)"}
NDS=${6:?"n_datasets required"}
MAXPAR=${7:-30}
DS_START=${8:-1}
DS_END=$((DS_START + NDS - 1))

case "$STRUCTURE" in
  linCmt|ode) ;;
  *) echo "ERROR: structure must be linCmt or ode, got '$STRUCTURE'" >&2; exit 1 ;;
esac

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
JOBNAME="bench_N${SAMPLE_N}_scn${SCN_PAD}_${STRUCTURE}_${EST}_${OPT}"

LOGDIR="$REPO_ROOT/logs/bench_N${SAMPLE_N}"
mkdir -p "$LOGDIR"

echo "Submitting: ${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}"
echo "  REPO_ROOT   = $REPO_ROOT"
echo "  SCRIPTS_DIR = $SCRIPTS_DIR"
echo "  logs        = $LOGDIR/${JOBNAME}.<jobid>.<idx>.{out,err}"

# --- optional knobs (set via env) ------------------------------------------
#   DEPEND_JOBID : if set, this array waits until that job ENDS (done OR exit)
#                  before starting -> serialize arrays to avoid host contention.
#   EXCLUSIVE=1  : request the whole node (-x) so the 3 SCM fork workers are
#                  never starved for cores -> trustworthy wall-time measurement.
#   JOBID_FILE   : if set, the submitted array's job id is written here so a
#                  caller (submit_all_arrays.sh CHAIN mode) can chain on it.
EXTRA_BSUB=()
if [ -n "${DEPEND_JOBID:-}" ]; then
  EXTRA_BSUB+=( -w "ended(${DEPEND_JOBID})" )
  echo "  depends on  = ended(${DEPEND_JOBID})"
fi
if [ "${EXCLUSIVE:-0}" = "1" ]; then
  EXTRA_BSUB+=( -x )
  echo "  exclusive   = -x (whole node)"
fi

# capture bsub stdout so we can extract the job id
bsub_out="$(
  bsub \
    "${EXTRA_BSUB[@]}" \
    -J "${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}" \
    -o "${LOGDIR}/${JOBNAME}.%J.%I.out" \
    -e "${LOGDIR}/${JOBNAME}.%J.%I.err" \
    -env "all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, SAMPLE_N=${SAMPLE_N}, SCN=${SCN}, EST=${EST}, OPT=${OPT}, STRUCTURE=${STRUCTURE}" \
    < "$HERE/bench_array.lsf"
)"
echo "$bsub_out"

# bsub prints:  Job <261262> is submitted to queue <short>.
JOBID="$(printf '%s\n' "$bsub_out" | sed -n 's/^Job <\([0-9]\+\)>.*/\1/p' | head -n1)"
if [ -n "${JOBID_FILE:-}" ] && [ -n "$JOBID" ]; then
  printf '%s\n' "$JOBID" > "$JOBID_FILE"
fi
echo "Submitted job id: ${JOBID:-<unparsed>}"
