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
# Optional: also stage each parsed record into the aggregator layout
#   ${STAGE_DST}/N<N>/scn<SS>_<STRUCT>/<EST>_<OPT>/res_ds<D>.rds
# so no separate restage loop is needed.  Leave empty to only write records/.
#   e.g.  STAGE_DST=output/psn_scm_full0727/ResforAggregation1 FORCE=1 bash parse_full_grid.sh
STAGE_DST="${STAGE_DST:-}"
STAGE_STRUCT="${STAGE_STRUCT:-advan4}"
STAGE_EST="${STAGE_EST:-nonmem_scm}"
STAGE_OPT="${STAGE_OPT:-focei}"
# STAGE_ONLY=1 writes res_ds<D>.rds DIRECTLY to STAGE_DST and skips records/
# entirely (parser --out_rds points at the res path).  Requires STAGE_DST.
STAGE_ONLY="${STAGE_ONLY:-0}"
# STRUCT filter: when set (e.g. STRUCT=ode), parse/stage ONLY manifest rows whose
# structure column matches -- so an ODE sweep can be targeted without touching
# advan4 rows in a shared manifest.  Empty = all structures.  The manifest's
# per-row structure (7th column) overrides STAGE_STRUCT for the res-layout path,
# so advan4 and ode land in scnSS_advan4/ vs scnSS_ode/ correctly.
STRUCT="${STRUCT:-}"

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
# manifest columns: N,scenario,dataset,cell,jobid,submitted_utc[,structure]
# structure is optional (older manifests lack it) -> default advan4.
while IFS=, read -r N SCEN ds cell jobid ts struct; do
  [ "${N}" = "N" ] && continue                       # header
  [ -n "${cell:-}" ] || continue
  struct="${struct:-advan4}"                          # back-compat: no col -> advan4
  # STRUCT filter: skip rows of other structures when targeting one (e.g. ode)
  [ -n "${STRUCT}" ] && [ "${struct}" != "${STRUCT}" ] && continue
  # per-cell res-layout structure segment (manifest struct wins over env default)
  cell_struct="${struct:-${STAGE_STRUCT}}"
  # record path mirrors the cell with runs/ -> records/ (same rule the parser uses)
  rec="$(printf '%s' "${cell}" | sed 's#/runs/#/records/#')/psn_scm_record.rds"

  # target res_ds<D>.rds path in the aggregator layout (used by both STAGE modes)
  if [ -n "${STAGE_DST}" ]; then
    scn2=$(printf '%02d' "${SCEN}")                  # zero-pad scenario -> scnSS
    dsn=$((10#${ds}))                                # strip any zero-pad -> res_ds<int>
    stage_dir="${STAGE_DST}/N${N}/scn${scn2}_${cell_struct}/${STAGE_EST}_${STAGE_OPT}"
    res="${stage_dir}/res_ds${dsn}.rds"
  fi

  if [ ! -f "${cell}/logs/timing.json" ]; then
    n_wait=$((n_wait + 1)); continue                 # not finished yet
  fi
  # existence check: STAGE_ONLY checks the res file (no records written); else record
  if [ "${FORCE}" != "1" ]; then
    if [ "${STAGE_ONLY}" = "1" ] && [ -n "${STAGE_DST}" ]; then
      [ -f "${res}" ] && { n_skip=$((n_skip + 1)); continue; }
    elif [ -f "${rec}" ]; then
      n_skip=$((n_skip + 1)); continue               # already parsed
    fi
  fi

  # STAGE_ONLY: write res_ds<D>.rds DIRECTLY (no records/ copy).  Else write the
  # record under records/ (default) and optionally copy it to the res layout.
  if [ "${STAGE_ONLY}" = "1" ] && [ -n "${STAGE_DST}" ]; then
    mkdir -p "${stage_dir}"
    parse_out=(--out_rds "${res}")
  else
    parse_out=()
  fi

  if "${RSCRIPT_BIN}" script/psn_scm/parse_psn_scm.R \
        --cell "${cell}" --N "${N}" --scenario "${SCEN}" --dataset "${ds}" \
        --structure "${cell_struct}" \
        "${parse_out[@]}" \
        >/dev/null 2>&1; then
    n_ok=$((n_ok + 1))
    # STAGE_DST (copy mode): stage the fresh record into the aggregator layout
    if [ "${STAGE_ONLY}" != "1" ] && [ -n "${STAGE_DST}" ] && [ -f "${rec}" ]; then
      mkdir -p "${stage_dir}"
      cp "${rec}" "${res}"
    fi
  else
    echo "  WARN: parse failed for ${cell}" >&2
    n_fail=$((n_fail + 1))
  fi
done < "${MANIFEST}"

echo "parsed: ${n_ok}  skipped(existing): ${n_skip}  waiting(unfinished): ${n_wait}  failed: ${n_fail}"
# report the structure segment actually used (filter value, else per-cell/default)
struct_shown="${STRUCT:-${STAGE_STRUCT}}"
if [ "${STAGE_ONLY}" = "1" ] && [ -n "${STAGE_DST}" ]; then
  echo "res     -> ${STAGE_DST}/N*/scn*_${struct_shown}/${STAGE_EST}_${STAGE_OPT}/res_ds*.rds  (records/ skipped)"
else
  echo "records -> ${BENCH_ROOT}/records/"
  [ -n "${STAGE_DST}" ] && echo "staged  -> ${STAGE_DST}/N*/scn*_${struct_shown}/${STAGE_EST}_${STAGE_OPT}/res_ds*.rds"
fi
