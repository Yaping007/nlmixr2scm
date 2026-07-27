#!/usr/bin/env bash
# ==============================================================================
# refit_winner.sh  --  standalone STAGE-2 covariance refit for a finished cell
# ------------------------------------------------------------------------------
# Salvages cells whose `scm` search already completed but whose refit failed
# (e.g. PsN 5.5 pharmpy-missing bug that corrupts scm_dir/final_models/).  Reuses
# the completed scm_dir + logs/scmlog.txt; does NOT re-run the search.
#
# Run from the REPO ROOT, pass one or more cell dirs:
#   bash script/psn_scm/refit_winner.sh output/psn_scm/runs/N300/scn16/ds00{1..5}
#
# For each cell it picks the real winning candidate model (matched to the final
# selected set in scmlog.txt), appends $COVARIANCE, executes it locally, and
# writes logs/refit.{lst,ext} + updates logs/timing.json (adds refit_sec/cpu and
# flips scope->scm+refit, refit_ok->1).  Best-effort: warns and skips on failure.
# ==============================================================================
set -uo pipefail   # NOT -e: a single cell failure must not abort the batch

NM_VERSION="${NM_VERSION:-7.5.1-gfortran}"
PSN_MODULE="${PSN_MODULE:-PsN/5.5.0-GCCcore-12.3.0}"
RUN_ON_LSF_FLAG="${RUN_ON_LSF_FLAG:--no-run_on_lsf}"
REFIT_COV_FILE_NAME="${REFIT_COV_FILE_NAME:-final_refit_cov.txt}"

# absolute path to seed_refit.R (refit_one runs inside `cd "${cell}"`).
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SEED_REFIT="${SEED_REFIT:-${REPO_ROOT}/script/psn_scm/seed_refit.R}"

export TMPDIR="${TMPDIR:-/scratch/$(whoami)}"; mkdir -p "${TMPDIR}" 2>/dev/null || true

# ---- module load (same self-contained pattern as submit_scm.sh) ------------
if ! type module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
              /usr/share/lmod/lmod/init/bash "${MODULESHOME:-/usr/share/Modules}/init/bash"; do
    [ -f "${init}" ] && { . "${init}"; break; }
  done
fi
module purge; module load "${PSN_MODULE}"
command -v execute >/dev/null 2>&1 || { echo "ERROR: execute not on PATH" >&2; exit 127; }

# seed_refit.R needs an R interpreter, but the R module (gomkl toolchain) differs
# from PsN's (GCCcore-12.3.0); loading it in THIS shell would swap dependencies
# and could break execute.  So run seed_refit in an ISOLATED SUBSHELL that loads
# its own R module and leaves the parent PsN environment untouched.
# Override with RSCRIPT=/abs/path/Rscript (skips module load) or R_MODULE=<name>.
R_MODULE="${R_MODULE:-R/4.3.1-gomkl-2022a-0.1}"
run_seed() {  # args passed straight to seed_refit.R
  if [ -n "${RSCRIPT:-}" ]; then
    "${RSCRIPT}" "${SEED_REFIT:-${REPO_ROOT}/script/psn_scm/seed_refit.R}" "$@"; return $?
  fi
  ( module purge >/dev/null 2>&1 || true
    module load "${R_MODULE}" >/dev/null 2>&1 || true
    command -v Rscript >/dev/null 2>&1 || { echo "no Rscript in ${R_MODULE}" >&2; exit 127; }
    Rscript "${SEED_REFIT:-${REPO_ROOT}/script/psn_scm/seed_refit.R}" "$@" )
}
echo "seed_refit R module -> ${RSCRIPT:-${R_MODULE}}"

[ "$#" -ge 1 ] || { echo "usage: $0 <cell_dir> [cell_dir ...]" >&2; exit 2; }

# tag+state signature of a candidate .mod: for each ";;;TAG-DEFINITION START..END"
# block, classify the shape -> PsN state (power=5, exp=4, cat=2) and print
# "TAG-STATE", sorted.  This disambiguates models that share covariate NAMES but
# carry one at a different state (e.g. forward CLBW at exp/4 vs winner power/5).
mod_sig() {
  awk '
    /-DEFINITION START/ { tag=$0; sub(/.*;;;[[:space:]]*/,"",tag);
                          sub(/-DEFINITION START.*/,"",tag); intag=1; body=""; next }
    /-DEFINITION END/   { s=0;
        if (body ~ /\*\*THETA/)                     s=5;
        else if (body ~ /EXP\(THETA/)               s=4;
        else if (body ~ /1[[:space:]]*\+[[:space:]]*THETA/) s=2;
        print tag"-"s; intag=0; next }
    intag { body=body" "$0 }
  ' "$1" | sort | tr '\n' ' ' | sed 's/ *$//'
}

refit_one() {
  local cell="$1"
  [ -d "${cell}/scm_dir" ] || { echo "SKIP ${cell}: no scm_dir"; return; }
  echo "=== refit ${cell} ==="
  ( cd "${cell}" || return
    # WIPE any stale refit artifacts first.  PsN `execute -directory=refit_run
    # -clean=0` will REUSE an existing refit_run/NM_run1 (from an earlier base-
    # model attempt) and copy its OLD results back WITHOUT re-running -- which
    # silently yields a base-model .lst/.ext (5 thetas, wrong OFV).  A fresh dir
    # forces a real compile+run of the seeded winner.
    rm -rf refit
    mkdir -p logs refit
    [ -f logs/scmlog.txt ] || cp -f scm_dir/scmlog.txt logs/scmlog.txt 2>/dev/null || true

    # final selected set WITH state: last forward state per (param,covar), minus
    # any relation moved to state 1 (removed) in a backward step.
    kept_sig="$(awk '
      /relation chosen in this forward step:/  {n=split($NF,a,"-"); st[a[1] a[2]]=a[n]}
      /relation chosen in this backward step:/ {n=split($NF,a,"-");
                                                if(a[n]==1) delete st[a[1] a[2]];
                                                else        st[a[1] a[2]]=a[n]}
      END{ for(t in st) print t"-"st[t] }
    ' logs/scmlog.txt 2>/dev/null | sort | tr '\n' ' ' | sed 's/ *$//' || true)"
    echo "  final selected (tag-state): ${kept_sig}"

    local winner="" m first sig
    while IFS= read -r m; do
      [ -f "${m}" ] || continue
      first="$(grep -vE '^\s*(;|$)' "${m}" | head -n1 || true)"
      case "${first}" in \$PROB*|\$SIZES*) : ;; *) continue ;; esac
      sig="$(mod_sig "${m}")"
      [ "${sig}" = "${kept_sig}" ] && { winner="${m}"; break; }
    done < <(find scm_dir -name '*.mod' -not -path '*/NM_run*' \
                 -not -path '*/temp_dir*' -not -path '*/final_models/*' 2>/dev/null)

    if [ -z "${winner}" ]; then
      echo "  WARN: no exact tag-state match; refusing base fallback -> skip"
      return
    fi
    echo "  winner: ${winner}"

    # Seed the refit from the winner's .ext final estimates (update_inits
    # replacement for the missing pharmpy) and force MAXEVAL=0, so NONMEM
    # computes the covariance AT the winner and the refit OFV is IDENTICAL to the
    # SCM winner OFV.  Without this, execute re-minimises from screening initials
    # (covariate thetas ~0.001) and drifts to a different, wrong optimum.
    local winner_ext="${winner%.mod}.ext"
    if [ ! -f "${winner_ext}" ]; then
      echo "  WARN: winner .ext missing (${winner_ext}); skip"; return
    fi
    run_seed --winner_mod "${winner}" --winner_ext "${winner_ext}" \
        --out refit/refit.mod --data "../data.csv" \
        --cov_file "${REFIT_COV_FILE_NAME}" --maxeval 0 \
      || { echo "  WARN: seed_refit failed; skip"; return; }

    if ( cd refit && /usr/bin/time -v -o ../logs/time_refit_raw.txt \
           execute refit.mod -directory=refit_run -nm_version="${NM_VERSION}" \
                   ${RUN_ON_LSF_FLAG} -clean=0 ) 2>&1 | tee logs/refit_console.log
    then
      cp -f refit/refit.lst logs/refit.lst 2>/dev/null || cp -f refit/refit_run/refit.lst logs/refit.lst 2>/dev/null || true
      cp -f refit/refit.ext logs/refit.ext 2>/dev/null || cp -f refit/refit_run/refit.ext logs/refit.ext 2>/dev/null || true
      # sanity: with MAXEVAL=0 the refit OFV must equal the SCM winner OFV.
      wofv="$(awk '$1==-1000000000{print $NF}' "${winner_ext}" 2>/dev/null | tail -1)"
      rofv="$(awk '$1==-1000000000{print $NF}' logs/refit.ext 2>/dev/null | tail -1)"
      if [ -n "${wofv}" ] && [ -n "${rofv}" ]; then
        d="$(awk -v a="${wofv}" -v b="${rofv}" 'BEGIN{d=a-b; if(d<0)d=-d; print d}')"
        awk -v d="${d}" 'BEGIN{exit !(d>0.01)}' \
          && echo "  WARN: refit OFV ${rofv} != winner OFV ${wofv} (|d|=${d}) -- seeding/exec suspect" \
          || echo "  OFV check OK: refit ${rofv} == winner ${wofv}"
      fi
      rcpu="$(awk -F': ' '/User time/{u=$2}/System time/{s=$2}END{printf "%.3f",u+s}' logs/time_refit_raw.txt)"
      # patch timing.json: scope->scm+refit, refit_ok->1, add refit_cpu_sec
      if [ -f logs/timing.json ]; then
        sed -i -E 's/"scope": *"[^"]*"/"scope": "scm+refit"/; s/"refit_ok": *[0-9]+/"refit_ok": 1/' logs/timing.json
        grep -q '"refit_cpu_sec"' logs/timing.json || \
          sed -i "s/\"cpu_sec\": \\([0-9.]*\\),/\"cpu_sec\": \\1,\n  \"refit_cpu_sec\": ${rcpu},/" logs/timing.json
      fi
      echo "  refit OK (cpu=${rcpu}s)"
    else
      echo "  WARN: refit execute failed"
    fi
  )
}

for cell in "$@"; do refit_one "${cell}"; done
echo "done."
