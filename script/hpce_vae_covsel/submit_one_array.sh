#!/bin/bash
# ============================================================================
# submit_one_array.sh
# ----------------------------------------------------------------------------
# Submit ONE LSF job array for one (N, scenario, structure) cell of the VAE
# covariate-selection pilot, with array indices 1..n_datasets (index = dataset).
#
# Usage:
#   bash script/hpce_vae_covsel/submit_one_array.sh <N> <scn> <structure> <n_ds> [maxpar] [ds_start]
#
# Example (pilot probe: N80, scn16, linCmt, 5 datasets, 20 parallel):
#   bash script/hpce_vae_covsel/submit_one_array.sh 80 16 linCmt 5 20
# ============================================================================

set -euo pipefail

SAMPLE_N=${1:?"Usage: $0 <N:40|80|300> <scn:1..16> <structure:linCmt|ode> <n_ds> [maxpar] [ds_start]"}
SCN=${2:?"scenario required"}
STRUCTURE=${3:?"structure required (linCmt|ode)"}
NDS=${4:?"n_datasets required"}
MAXPAR=${5:-20}
DS_START=${6:-1}
DS_END=$((DS_START + NDS - 1))

case "$STRUCTURE" in
  linCmt|ode) ;;
  *) echo "ERROR: structure must be linCmt or ode (got '$STRUCTURE')" >&2; exit 1 ;;
esac

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$HERE/../.." && pwd)}"

# Auto-detect scripts subdir holding the driver
SCRIPTS_DIR=""
for c in script/hpce_vae_covsel script scripts R; do
  if [ -f "$REPO_ROOT/$c/vae_covsel_driver.R" ]; then
    SCRIPTS_DIR="$c"; break
  fi
done
if [ -z "$SCRIPTS_DIR" ]; then
  echo "ERROR: vae_covsel_driver.R not found under $REPO_ROOT" >&2
  exit 1
fi

if [ ! -f "$HERE/vae_covsel_array.lsf" ]; then
  echo "ERROR: $HERE/vae_covsel_array.lsf missing" >&2
  exit 1
fi

OUT_ROOT="${OUT_ROOT:-output/vae_covsel_pilot}"
SCN_PAD=$(printf '%02d' "$SCN")
JOBNAME="vaecov_N${SAMPLE_N}_scn${SCN_PAD}_${STRUCTURE}"

LOGDIR="$REPO_ROOT/logs/vaecov_N${SAMPLE_N}"
mkdir -p "$LOGDIR"

echo "Submitting: ${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}"
echo "  REPO_ROOT   = $REPO_ROOT"
echo "  SCRIPTS_DIR = $SCRIPTS_DIR"
echo "  OUT_ROOT    = $OUT_ROOT"
echo "  logs        = $LOGDIR/${JOBNAME}.<jobid>.<idx>.{out,err}"

bsub \
  -J "${JOBNAME}[${DS_START}-${DS_END}]%${MAXPAR}" \
  -o "${LOGDIR}/${JOBNAME}.%J.%I.out" \
  -e "${LOGDIR}/${JOBNAME}.%J.%I.err" \
  -env "all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, SAMPLE_N=${SAMPLE_N}, SCN=${SCN}, STRUCTURE=${STRUCTURE}, OUT_ROOT=${OUT_ROOT}" \
  < "$HERE/vae_covsel_array.lsf"
