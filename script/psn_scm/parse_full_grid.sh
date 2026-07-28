#!/usr/bin/env bash
# ==============================================================================
# parse_full_grid.sh  --  parse every finished cell of the full grid -> records/
# ------------------------------------------------------------------------------
# Walks the manifest (written by run_full_grid.sh) and runs parse_psn_scm.R on
# each finished cell (has logs/timing.json), writing the light record next to it
# under records/ (runs/ -> records/).  Defaults to the isolated full-launch tree
# output/psn_scm_full0727/ (override with BENCH_ROOT= or MANIFEST=).
#
# Run from the REPO ROOT after `bjobs` is empty:
#
#     bash script/psn_scm/parse_full_grid.sh
#     FORCE=1 bash script/psn_scm/parse_full_grid.sh     # re-parse even if record exists
#
# Env overrides: RSCRIPT, R_MODULE (R for parsing), FORCE, MANIFEST.
# parse_psn_scm.R needs Rscript (not on PsN PATH); resolved once here.
# ==============================================================================
set -uo pipefail

R_MODULE="${R_MODULE:-R/4.3.1-gomkl-2022a-0.1}"
FORCE="${FORCE:-0}"
BENCH_ROOT="${BENCH_ROOT:-output/psn_scm_full0727}"
MANIFEST="${MANIFEST:-${BENCH_ROOT}/manifest.csv}"

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${REPO_ROOT}"

[ -f "${MANIFEST}" ] || { echo "ERROR: no manifest at ${MANIFEST}" >&2; exit 2; }

# ---- resolve Rscript for parsing --------------------------------------------
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
[ -n "${RSCRIPT_BIN}" ] || { echo "ERROR: no Rscript (set RSCRIPT / R_MODULE)." >&2; exit 127; }
echo "parse Rscript -> ${RSCRIPT_BIN}"

n_ok=0; n_skip=0; n_wait=0; n_fail=0
# manifest columns: N,scenario,dataset,cell,jobid,submitted_utc
while IFS=, read -r N SCEN ds cell jobid ts; do
  [ "${N}" = "N" ] && continue                       # header
  [ -n "${cell:-}" ] || continue
  # record path mirrors the cell with runs/ -> records/ (same rule the parser uses)
  rec="$(printf '%s' "${cell}" | sed 's#/runs/#/records/#')/psn_scm_record.rds"

  if [ ! -f "${cell}/logs/timing.json" ]; then
    n_wait=$((n_wait + 1)); continue                 # not finished yet
  fi
  if [ "${FORCE}" != "1" ] && [ -f "${rec}" ]; then
    n_skip=$((n_skip + 1)); continue                 # already parsed
  fi

  if "${RSCRIPT_BIN}" script/psn_scm/parse_psn_scm.R \
        --cell "${cell}" --N "${N}" --scenario "${SCEN}" --dataset "${ds}" \
        >/dev/null 2>&1; then
    n_ok=$((n_ok + 1))
  else
    echo "  WARN: parse failed for ${cell}" >&2
    n_fail=$((n_fail + 1))
  fi
done < "${MANIFEST}"

echo "parsed: ${n_ok}  skipped(existing): ${n_skip}  waiting(unfinished): ${n_wait}  failed: ${n_fail}"
echo "records -> output/psn_scm/records/"
