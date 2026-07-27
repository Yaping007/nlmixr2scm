#!/usr/bin/env bash
# ==============================================================================
# run_pilot_scn16_N300.sh  --  pilot sweep: scn16, N=300, datasets 1..5, w/ refit
# ------------------------------------------------------------------------------
# Exports 5 cells and submits each as its own single-node LSF job (base+search
# +mandatory covariance refit).  Run from the REPO ROOT on DaVinci HPCE.
#
#   bash script/psn_scm/run_pilot_scn16_N300.sh
#
# Each cell -> output/psn_scm/runs/N300/scn16/ds00{1..5}/.  After the jobs
# finish, parse every cell into output/psn_scm/records/ (see the parse loop at
# the bottom, or run it separately once bjobs is empty).
#
# Env overrides: N, SCEN, DS_LIST, CORES, WALL_MIN, RSCRIPT.
# ==============================================================================
set -euo pipefail

N="${N:-300}"
SCEN="${SCEN:-16}"
DS_LIST="${DS_LIST:-1 2 3 4 5}"
CORES="${CORES:-4}"
WALL_MIN="${WALL_MIN:-10000}"
RSCRIPT="${RSCRIPT:-Rscript}"

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${REPO_ROOT}"
RUNS_ROOT="output/psn_scm/runs"

scn2="$(printf '%02d' "${SCEN}")"

echo "=== PsN SCM pilot: N=${N} scn=${scn2} datasets=[${DS_LIST}] (scm+refit) ==="

for ds in ${DS_LIST}; do
  ds3="$(printf '%03d' "${ds}")"
  cell="${RUNS_ROOT}/N${N}/scn${scn2}/ds${ds3}"

  # 1) export inputs (data.csv, base.mod, run.scm, final_refit_cov.txt, submit_scm.sh)
  echo "--- export ${cell} ---"
  "${RSCRIPT}" script/psn_scm/export_one_dataset.R \
      --N "${N}" --scenario "${SCEN}" --dataset "${ds}"

  # 2) normalise line endings on the shell script (repo may carry CRLF)
  sed -i 's/\r$//' "${cell}/submit_scm.sh" 2>/dev/null || true

  # 3) submit one single-node job per cell
  ( cd "${cell}" && \
    bsub -n "${CORES}" -W "${WALL_MIN}" \
         -o "lsf.%J.out" -e "lsf.%J.err" \
         -J "psnscm_N${N}_s${scn2}_d${ds3}" \
         "bash submit_scm.sh" )
done

echo
echo "Submitted ${DS_LIST// /,} datasets.  Watch with:  bjobs -w"
echo
echo "When all jobs finish, parse each cell into records/:"
echo
for ds in ${DS_LIST}; do
  ds3="$(printf '%03d' "${ds}")"
  echo "  ${RSCRIPT} script/psn_scm/parse_psn_scm.R \\"
  echo "      --cell ${RUNS_ROOT}/N${N}/scn${scn2}/ds${ds3} \\"
  echo "      --N ${N} --scenario ${SCEN} --dataset ${ds}"
done
