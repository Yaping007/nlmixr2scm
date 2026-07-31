#!/usr/bin/env bash
# ==============================================================================
# run_full_grid.sh  --  full PsN SCM benchmark sweep on DaVinci HPCE
# ------------------------------------------------------------------------------
# Exports and submits one single-node LSF job per cell over the FULL grid:
#
#     N        in {40, 80, 300}          (cohorts)
#     scenario in {1 .. 16}              (all scenarios)
#     dataset  in {1 .. 100}             (ds001 .. ds100)
#
#   => 3 x 16 x 100 = 4800 cells, each running scm(base+search) + covariance
#      refit (submit_scm.sh) and writing its own logs/timing.json.
#
# Output goes to an ISOLATED tree, output/psn_scm_full0727/{runs,records}, kept
# SEPARATE from the working pilot (output/psn_scm/) so the full sweep never
# overwrites pilot cells.  Override with BENCH_ROOT=...
#
# Run from the REPO ROOT on DaVinci HPCE:
#
#     bash script/psn_scm/run_full_grid.sh                 # submit everything
#     N_LIST="300" SCEN_LIST="16" DS_MAX=10 \              # smoke a slice
#         bash script/psn_scm/run_full_grid.sh
#     DRY_RUN=1 bash script/psn_scm/run_full_grid.sh       # print, don't submit
#
# Env overrides:
#   N_LIST     cohorts               (default "40 80 300")
#   SCEN_LIST  scenarios             (default "1 2 ... 16")
#   DS_MIN/DS_MAX  dataset range     (default 1 / 100)
#   CORES      -n per job            (default 4)
#   WALL_MIN   -W minutes per job    (default 10000)
#   QUEUE      -q LSF queue          (default: unset -> site default)
#   THROTTLE   max pending+running jobs for this user; the driver paces
#              submission to stay under it                (default 400)
#   RSCRIPT    Rscript path; else R_MODULE is module-loaded for the export step
#   R_MODULE   R module for export   (default R/4.3.1-gomkl-2022a-0.1)
#   BENCH_ROOT output tree           (default output/psn_scm_full0727)
#   DRY_RUN=1  print the plan, submit nothing
#
# NOTE: the export step needs Rscript, which is NOT on the PsN module PATH.  The
# driver resolves R once (RSCRIPT or R_MODULE) and uses it only for export;
# submit_scm.sh loads PsN itself inside each job.
#
# A manifest (one row per submitted cell: keys, jobid, cell path) is appended to
# output/psn_scm/manifest.csv for tracking ~4800 cells and locating failures.
# ==============================================================================
set -uo pipefail   # NOT -e: one failed export/bsub must not abort the sweep

N_LIST="${N_LIST:-40 80 300}"
SCEN_LIST="${SCEN_LIST:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16}"
DS_MIN="${DS_MIN:-1}"
DS_MAX="${DS_MAX:-100}"
CORES="${CORES:-4}"
WALL_MIN="${WALL_MIN:-10000}"
QUEUE="${QUEUE:-}"
THROTTLE="${THROTTLE:-400}"
R_MODULE="${R_MODULE:-R/4.3.1-gomkl-2022a-0.1}"
DRY_RUN="${DRY_RUN:-0}"
# STRUCT selects the structural model for the sweep (advan4|ode); run ode into a
# separate BENCH_ROOT.  STRUCT_TOL_ARGS forwards --screen_tol/--refit_tol etc.
STRUCT="${STRUCT:-advan4}"
STRUCT_TOL_ARGS="${STRUCT_TOL_ARGS:-}"

# Isolated output tree for the FULL launch -- kept SEPARATE from the pilot
# (output/psn_scm/runs, records) so a full sweep never overwrites the working
# pilot cells.  Override with BENCH_ROOT=... if needed.
BENCH_ROOT="${BENCH_ROOT:-output/psn_scm_full0727}"

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${REPO_ROOT}"
RUNS_ROOT="${BENCH_ROOT}/runs"
MANIFEST="${BENCH_ROOT}/manifest.csv"
mkdir -p "$(dirname "${MANIFEST}")"

# ---- resolve an Rscript for the export step (login shell, no PsN module) -----
if ! type module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
              /usr/share/lmod/lmod/init/bash "${MODULESHOME:-/usr/share/Modules}/init/bash"; do
    [ -f "${init}" ] && { . "${init}"; break; }
  done
fi
resolve_rscript() {
  if [ -n "${RSCRIPT:-}" ] && command -v "${RSCRIPT}" >/dev/null 2>&1; then
    echo "${RSCRIPT}"; return; fi
  if command -v Rscript >/dev/null 2>&1; then command -v Rscript; return; fi
  module load "${R_MODULE}" >/dev/null 2>&1 && command -v Rscript >/dev/null 2>&1 && {
    command -v Rscript; return; }
}
RSCRIPT_BIN="$(resolve_rscript || true)"
if [ -z "${RSCRIPT_BIN}" ]; then
  echo "ERROR: no Rscript for export (set RSCRIPT=/path or R_MODULE=<module>)." >&2
  exit 127
fi
echo "export Rscript -> ${RSCRIPT_BIN}"

# ---- throttle: pause while my pending+running jobs exceed THROTTLE -----------
# Count PEND/RUN jobs for this user.  Robust to bjobs printing nothing (no jobs):
# pipe through wc -l so the result is ALWAYS a single integer (grep -c on empty
# input exits non-zero and, combined with `|| echo 0`, used to emit "0\n0" which
# broke the numeric test below).  `bjobs -noheader` avoids counting the header.
my_jobs() {
  local n
  n="$(bjobs -noheader -o 'stat' 2>/dev/null | grep -cE '^(PEND|RUN)$')"
  # strip any stray whitespace/newlines; default to 0 if empty
  n="$(printf '%s' "${n}" | tr -dc '0-9')"
  printf '%s' "${n:-0}"
}
wait_for_slot() {
  while :; do
    local n; n="$(my_jobs)"
    [ "${n:-0}" -lt "${THROTTLE}" ] && break
    echo "  throttle: ${n} jobs in queue (>= ${THROTTLE}); sleeping 60s..."
    sleep 60
  done
}

# ---- manifest header (once) --------------------------------------------------
[ -f "${MANIFEST}" ] || echo "N,scenario,dataset,cell,jobid,submitted_utc,structure" > "${MANIFEST}"

n_cells=0; n_submitted=0; n_export_fail=0
echo "=== PsN SCM full grid: N=[${N_LIST}] scn=[${SCEN_LIST}] ds=${DS_MIN}..${DS_MAX} ==="
[ "${DRY_RUN}" = "1" ] && echo "(DRY_RUN: printing plan, submitting nothing)"

for N in ${N_LIST}; do
  for SCEN in ${SCEN_LIST}; do
    scn2="$(printf '%02d' "${SCEN}")"
    for ds in $(seq "${DS_MIN}" "${DS_MAX}"); do
      ds3="$(printf '%03d' "${ds}")"
      cell="${RUNS_ROOT}/N${N}/scn${scn2}/ds${ds3}"
      n_cells=$((n_cells + 1))

      if [ "${DRY_RUN}" = "1" ]; then
        echo "  would export+submit ${cell}"
        continue
      fi

      # 1) export inputs
      if ! "${RSCRIPT_BIN}" script/psn_scm/export_one_dataset.R \
             --N "${N}" --scenario "${SCEN}" --dataset "${ds}" \
             --structure "${STRUCT}" ${STRUCT_TOL_ARGS} \
             --out_root "${RUNS_ROOT}" >/dev/null 2>&1; then
        echo "  WARN: export failed for ${cell}; skipping" >&2
        n_export_fail=$((n_export_fail + 1))
        continue
      fi

      # 2) strip CRLF on the shell wrapper (repo may carry CRLF)
      sed -i 's/\r$//' "${cell}/submit_scm.sh" 2>/dev/null || true

      # 3) throttle, then submit one single-node job per cell
      wait_for_slot
      jobname="psnscm_N${N}_s${scn2}_d${ds3}"
      qflag=""; [ -n "${QUEUE}" ] && qflag="-q ${QUEUE}"
      out="$( cd "${cell}" && \
        bsub -n "${CORES}" -W "${WALL_MIN}" ${qflag} \
             -o "lsf.%J.out" -e "lsf.%J.err" -J "${jobname}" \
             "bash submit_scm.sh" 2>&1 )"
      echo "  ${cell}: ${out}"
      jobid="$(printf '%s' "${out}" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+' | head -1)"
      echo "${N},${SCEN},${ds},${cell},${jobid:-NA},$(date -u +%Y-%m-%dT%H:%M:%SZ),${STRUCT}" >> "${MANIFEST}"
      n_submitted=$((n_submitted + 1))
    done
  done
done

echo
echo "grid cells: ${n_cells}  submitted: ${n_submitted}  export-fail: ${n_export_fail}"
echo "manifest -> ${MANIFEST}"
echo "watch:  bjobs -w    |    summary:  bjobs -o 'stat' | sort | uniq -c"
echo
echo "When all jobs finish, parse the whole grid into records/:"
echo "  bash script/psn_scm/parse_full_grid.sh   # (loops manifest -> parse_psn_scm.R)"
