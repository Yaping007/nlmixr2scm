#!/usr/bin/env bash
# ==============================================================================
# run_full_grid_array.sh  --  full PsN SCM benchmark sweep as ONE LSF job array
# ------------------------------------------------------------------------------
# Same grid as run_full_grid.sh, but submits a SINGLE LSF *job array* instead of
# thousands of individual bsub calls:
#
#     N        in {40, 80, 300}          (cohorts)
#     scenario in {1 .. 16}              (all scenarios)
#     dataset  in {1 .. 100}             (ds001 .. ds100)
#
#   => up to 3 x 16 x 100 = 4800 cells = one array  psnscm[1-4800]%THROTTLE
#
# Why an array (vs run_full_grid.sh's per-cell bsub):
#   * ONE submission, not ~4800 -- no scheduler hammering, no bjobs clutter.
#   * LSF paces concurrency itself via the % throttle (psnscm[1-N]%400); the
#     old wait_for_slot sleep-loop is gone.
#   * bkill <arrayid> cancels the whole sweep; bkill "<arrayid>[37]" one cell.
#   * memory + walltime set explicitly -> no "3000 MB default" warning, and
#     jobs land in a normal queue (pilot: 5-8 min/cell, ~110 MB peak).
#
# TWO PHASES:
#   Phase 1 (login node): export every cell's inputs (fast, no NONMEM) and write
#            an index file  ${BENCH_ROOT}/array_index.tsv  mapping
#            LSF array index -> cell dir.  Also (re)writes array_task.sh.
#   Phase 2: submit ONE array; each element reads its line from array_index.tsv
#            (by $LSB_JOBINDEX), cd's into that cell, runs submit_scm.sh.
#
# Output tree is the ISOLATED full-launch tree (kept separate from the pilot):
#     ${BENCH_ROOT}/{runs,records}   default output/psn_scm_full0727
#
# Run from the REPO ROOT on DaVinci HPCE:
#     bash script/psn_scm/run_full_grid_array.sh                 # full grid
#     N_LIST="300" SCEN_LIST="16" DS_MAX=5 \                     # smoke a slice
#         bash script/psn_scm/run_full_grid_array.sh
#     DRY_RUN=1 bash script/psn_scm/run_full_grid_array.sh       # export+index only
#
# Clean up the earlier per-cell singles first (run_full_grid.sh) if still queued:
#     bkill -J 'psnscm_*'
#
# Env overrides:
#   N_LIST     cohorts               (default "40 80 300")
#   SCEN_LIST  scenarios             (default "1 2 ... 16")
#   DS_MIN/DS_MAX  dataset range     (default 1 / 100)
#   CORES      -n per array element  (default 4)
#   WALL_MIN   -W minutes per element(default 60   -- pilot worst case ~8 min)
#   MEM_MB     -M / rusage[mem] MB   (default 2000 -- pilot peak ~110 MB)
#   THROTTLE   max concurrent array elements  (default 400  -> psnscm[1-N]%400)
#   QUEUE      -q LSF queue          (default: unset -> site default)
#   RSCRIPT / R_MODULE  R for the export step (NOT on the PsN module PATH)
#   BENCH_ROOT output tree           (default output/psn_scm_full0727)
#   JOBNAME    array job name        (default psnscm_full)
#   SKIP_DONE=1 skip cells that already have logs/timing.json (resume a partial
#              sweep -- resubmit only missing/failed cells)
#   DRY_RUN=1  export + build index, submit nothing
# ==============================================================================
set -uo pipefail   # NOT -e: one failed export must not abort the sweep

N_LIST="${N_LIST:-40 80 300}"
SCEN_LIST="${SCEN_LIST:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16}"
DS_MIN="${DS_MIN:-1}"
DS_MAX="${DS_MAX:-100}"
CORES="${CORES:-4}"
# STRUCT selects the NONMEM structural model for the WHOLE sweep:
#   "advan4" (default) -- ADVAN4 analytic (current head-to-head vs runSCM linCmt)
#   "ode"              -- ADVAN13 general-ODE, tol-matched to runSCM structure=ode
# Run ODE into a SEPARATE BENCH_ROOT so the advan4 sweep is never touched, e.g.
#   STRUCT=ode BENCH_ROOT=output/psn_scm_ode0729 bash run_full_grid_array.sh
# Optional --screen_tol/--refit_tol/--screen_atol/--refit_atol are forwarded to
# export via the STRUCT_TOL_ARGS passthrough (default: TOL=6 ATOL=8 both tiers).
STRUCT="${STRUCT:-advan4}"
STRUCT_TOL_ARGS="${STRUCT_TOL_ARGS:-}"
# Continuous-covariate valid_states menu forwarded to the exporter. Default
# "1,4,5" = nlmixr2's {power, lin} competing shapes; set CONTINUOUS_STATES=1,5
# for the power-only covariate space (into a SEPARATE BENCH_ROOT).
CONTINUOUS_STATES="${CONTINUOUS_STATES:-1,4,5}"
WALL_MIN="${WALL_MIN:-60}"
MEM_MB="${MEM_MB:-2000}"
THROTTLE="${THROTTLE:-400}"
QUEUE="${QUEUE:-}"
R_MODULE="${R_MODULE:-R/4.3.1-gomkl-2022a-0.1}"
JOBNAME="${JOBNAME:-psnscm_full}"
DRY_RUN="${DRY_RUN:-0}"
# SKIP_DONE=1 -> Phase 1 skips any cell that already has logs/timing.json
# (a completed run).  Lets you resubmit ONLY the missing/failed cells after a
# partial/aborted sweep without re-running the ones that already finished.
SKIP_DONE="${SKIP_DONE:-0}"

# Isolated output tree for the FULL launch -- SEPARATE from the pilot
# (output/psn_scm/) so the sweep never overwrites working pilot cells.
BENCH_ROOT="${BENCH_ROOT:-output/psn_scm_full0727}"

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${REPO_ROOT}"
RUNS_ROOT="${BENCH_ROOT}/runs"
MANIFEST="${BENCH_ROOT}/manifest.csv"
# Namespace the index/task per JOBNAME so concurrent (or aborted) launches into
# the SAME BENCH_ROOT cannot truncate each other's live array index.  A running
# array element resolves its cell from the index at runtime, so a shared
# array_index.tsv is a footgun: a second launch's Phase-1 ': > INDEX' would
# clobber the first array's mapping mid-flight.
INDEX="${BENCH_ROOT}/array_index.${JOBNAME}.tsv"
TASK="${BENCH_ROOT}/array_task.${JOBNAME}.sh"
LOGDIR="${BENCH_ROOT}/lsf"
mkdir -p "${RUNS_ROOT}" "${LOGDIR}"

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

# ==============================================================================
# Phase 1: export every cell + build the array index (idx -> cell)
# ==============================================================================
: > "${INDEX}"                                        # truncate index
[ -f "${MANIFEST}" ] || echo "N,scenario,dataset,cell,jobid,submitted_utc,structure" > "${MANIFEST}"

idx=0; n_export_fail=0
echo "=== Phase 1: export + index  N=[${N_LIST}] scn=[${SCEN_LIST}] ds=${DS_MIN}..${DS_MAX} ==="
for N in ${N_LIST}; do
  for SCEN in ${SCEN_LIST}; do
    scn2="$(printf '%02d' "${SCEN}")"
    for ds in $(seq "${DS_MIN}" "${DS_MAX}"); do
      ds3="$(printf '%03d' "${ds}")"
      cell="${RUNS_ROOT}/N${N}/scn${scn2}/ds${ds3}"

      # SKIP_DONE: a cell with logs/timing.json already finished -> don't
      # re-export or resubmit it (resume a partial/aborted sweep cheaply).
      if [ "${SKIP_DONE}" = "1" ] && [ -f "${cell}/logs/timing.json" ]; then
        continue
      fi

      # export inputs (fast, no NONMEM)
      if ! "${RSCRIPT_BIN}" script/psn_scm/export_one_dataset.R \
             --N "${N}" --scenario "${SCEN}" --dataset "${ds}" \
             --structure "${STRUCT}" ${STRUCT_TOL_ARGS} \
             --continuous_states "${CONTINUOUS_STATES}" \
             --out_root "${RUNS_ROOT}" >/dev/null 2>&1; then
        echo "  WARN: export failed for ${cell}; skipping" >&2
        n_export_fail=$((n_export_fail + 1))
        continue
      fi
      sed -i 's/\r$//' "${cell}/submit_scm.sh" 2>/dev/null || true

      idx=$((idx + 1))
      # index line:  idx <TAB> N <TAB> scenario <TAB> dataset <TAB> cell
      printf '%d\t%s\t%s\t%s\t%s\n' "${idx}" "${N}" "${SCEN}" "${ds}" "${cell}" >> "${INDEX}"
    done
  done
done
NTASK=${idx}
echo "exported ${NTASK} cells  (export-fail: ${n_export_fail})"
echo "index -> ${INDEX}"
[ "${NTASK}" -gt 0 ] || { echo "ERROR: nothing to submit (0 cells)." >&2; exit 3; }

# ==============================================================================
# Phase 2a: write the per-element task script
# ==============================================================================
# Each array element runs this: look up its line in array_index.tsv by
# $LSB_JOBINDEX, cd into that cell, run submit_scm.sh (which loads PsN itself).
cat > "${TASK}" <<'TASK_EOF'
#!/usr/bin/env bash
set -uo pipefail
# $1 = index file (array_index.tsv);  element selected by $LSB_JOBINDEX
INDEX_FILE="$1"
i="${LSB_JOBINDEX:?LSB_JOBINDEX not set (must run as an LSF array element)}"
line="$(awk -F'\t' -v k="${i}" '$1==k {print; exit}' "${INDEX_FILE}")"
[ -n "${line}" ] || { echo "no index row for element ${i}" >&2; exit 4; }
cell="$(printf '%s' "${line}" | cut -f5)"
[ -d "${cell}" ] || { echo "cell dir missing: ${cell}" >&2; exit 5; }
cd "${cell}" || exit 6
echo "[array ${i}] $(date -u +%FT%TZ)  cell=${cell}"
exec bash submit_scm.sh
TASK_EOF
sed -i 's/\r$//' "${TASK}"
chmod +x "${TASK}"
echo "task  -> ${TASK}"

# ==============================================================================
# Phase 2b: submit ONE array  psnscm_full[1-NTASK]%THROTTLE
# ==============================================================================
qflag=""; [ -n "${QUEUE}" ] && qflag="-q ${QUEUE}"
arr="${JOBNAME}[1-${NTASK}]%${THROTTLE}"

echo
echo "=== Phase 2: submit array  ${arr} ==="
echo "  -n ${CORES}  -W ${WALL_MIN}min  -M ${MEM_MB}MB  rusage[mem=${MEM_MB}]  ${qflag:-(default queue)}"

if [ "${DRY_RUN}" = "1" ]; then
  echo "(DRY_RUN: exported ${NTASK} cells + wrote index/task; NOT submitting.)"
  echo "to submit for real, re-run without DRY_RUN=1."
  exit 0
fi

out="$(
  bsub -J "${arr}" \
       -n "${CORES}" -W "${WALL_MIN}" \
       -M "${MEM_MB}" -R "rusage[mem=${MEM_MB}]" \
       ${qflag} \
       -o "${LOGDIR}/%J.%I.out" -e "${LOGDIR}/%J.%I.err" \
       "bash ${TASK} ${INDEX}" 2>&1
)"
echo "  ${out}"
arrayid="$(printf '%s' "${out}" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+' | head -1)"

# ---- append manifest rows (one per element) so parse_full_grid.sh works as-is
# manifest columns: N,scenario,dataset,cell,jobid,submitted_utc,structure
# jobid recorded as <arrayid>[<idx>] so a failed element is locatable.  The
# trailing structure column tells parse_full_grid.sh which model each cell used
# (advan4|ode) so parsing/staging tag it correctly.
now="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
while IFS=$'\t' read -r i N SCEN ds cell; do
  echo "${N},${SCEN},${ds},${cell},${arrayid:-NA}[${i}],${now},${STRUCT}" >> "${MANIFEST}"
done < "${INDEX}"

echo
echo "array job: ${arrayid:-NA}   elements: ${NTASK}   throttle: %${THROTTLE}"
echo "manifest -> ${MANIFEST}"
echo "watch:   bjobs -A ${arrayid:-<id>}            # array summary"
echo "         bjobs -a -J '${JOBNAME}' | head       # element states"
echo "cancel:  bkill ${arrayid:-<id>}               # whole sweep"
echo "         bkill '${arrayid:-<id>}[37]'          # one element"
echo
echo "When the array finishes, parse the whole grid into records/:"
echo "  bash script/psn_scm/parse_full_grid.sh"
