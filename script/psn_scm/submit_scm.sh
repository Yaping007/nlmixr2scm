#!/usr/bin/env bash
# ==============================================================================
# submit_scm.sh  --  timed PsN SCM pipeline (Option A: single-node benchmark)
# ------------------------------------------------------------------------------
# Reproduces the nlmixr2 SCM benchmark's TOTAL wall time (base fit + SCM search
# + final covariance refit).  PsN's `scm` fits the base model as its own node-0
# and then runs forward+backward selection in ONE call, so STAGE 1 already
# covers nlmixr2's t_base + t_scm.  STAGE 2 is the single tight-tol covariance
# refit of the winning model (nlmixr2's t_refit).  Both stages are timed as one
# block; the refit is MANDATORY (SE + condition number come from it).
#
#     nlmixr2:  t_total = t_base + t_scm + t_refit
#     PsN     :  t_total = t_scm(base+search) + t_refit
#
# Runs on a SINGLE multi-core node with candidates LOCAL (-no-run_on_lsf) so
# /usr/bin/time reaps every NONMEM child's CPU (accurate cpu_sec + hog),
# mirroring runSCM(workers=3, rxThreads=1).
#
# Usage (submit the whole pipeline as one multi-core LSF job):
#     module load PsN/5.5.0-GCCcore-12.3.0
#     bsub -n 4 -W 10000 "bash submit_scm.sh"
#
# Writes logs/timing.json: total wall_sec + cpu_sec + hog for the pipeline.
# ==============================================================================
set -euo pipefail

NM_VERSION="${NM_VERSION:-7.5.1-gfortran}"
THREADS="${THREADS:-3}"             # parallel NONMEM candidate jobs == runSCM workers=3
SCM_CFG="${SCM_CFG:-run.scm}"
PSN_MODULE="${PSN_MODULE:-PsN/5.5.0-GCCcore-12.3.0}"
RUN_ON_LSF_FLAG="${RUN_ON_LSF_FLAG:--no-run_on_lsf}"
REFIT_COV_FILE="${REFIT_COV_FILE:-final_refit_cov.txt}"

# seed_refit.R = update_inits replacement (seeds refit initials from winner .ext)
# submit_scm.sh runs from the cell dir, where export_one_dataset.R ships a copy.
if [ -n "${SEED_REFIT:-}" ]; then :
elif [ -f "./seed_refit.R" ]; then SEED_REFIT="./seed_refit.R"
else SEED_REFIT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/seed_refit.R"
fi

mkdir -p logs

# PsN/MPI gotcha (site doc): missing /machines hostfile -> set a private TMPDIR.
export TMPDIR="${TMPDIR:-/scratch/$(whoami)}"
mkdir -p "${TMPDIR}" 2>/dev/null || true

# ---- self-contained module load -------------------------------------------
# purge first: the login/compute node defaults to GCCcore/14.2.0, which
# conflicts with PsN 5.5's GCCcore/12.3.0.  purge + load pulls the matching
# GCCcore 12.3.0 (and NONMEM 7.5.1-gfortran shares the same toolchain).
#
# `module` is a SHELL FUNCTION, not an executable, and it is NOT defined in a
# non-interactive batch shell -> we must source the init script first, else the
# load is silently skipped and `scm` is not on PATH (exit 127).
if ! type module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
              /usr/share/lmod/lmod/init/bash "${MODULESHOME:-/usr/share/Modules}/init/bash"; do
    if [ -f "${init}" ]; then . "${init}"; break; fi
  done
fi
module purge
module load "${PSN_MODULE}"

# hard fail early with a clear message if the PsN tools are not resolvable
for tool in execute scm; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "ERROR: '${tool}' not on PATH after loading ${PSN_MODULE}." >&2
    echo "       module list:" >&2; module list 2>&1 || true
    exit 127
  fi
done
echo "scm -> $(command -v scm)   execute -> $(command -v execute)"

# seed_refit.R needs an R interpreter, but the R module (gomkl toolchain) differs
# from PsN's (GCCcore-12.3.0); loading it in THIS shell would swap dependencies
# and could break execute/scm.  So run seed_refit in an ISOLATED SUBSHELL that
# loads its own R module and leaves the parent PsN environment untouched.
# Override with RSCRIPT=/abs/path/Rscript (skips module load) or R_MODULE=<name>.
R_MODULE="${R_MODULE:-R/4.3.1-gomkl-2022a-0.1}"
run_seed() {  # args passed straight to seed_refit.R
  if [ -n "${RSCRIPT:-}" ]; then
    "${RSCRIPT}" "${SEED_REFIT}" "$@"; return $?
  fi
  ( module purge >/dev/null 2>&1 || true
    module load "${R_MODULE}" >/dev/null 2>&1 || true
    command -v Rscript >/dev/null 2>&1 || { echo "no Rscript in ${R_MODULE}" >&2; exit 127; }
    Rscript "${SEED_REFIT}" "$@" )
}

# ============================================================================
# Timed pipeline: scm(base+search) -> refit, measured as ONE block.
#   wall via date (whole pipeline); cpu by summing per-command /usr/bin/time
#   raw files (each NONMEM child is local, so its CPU is captured).
# ============================================================================
start_iso="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
wall_start="$(date +%s.%N)"

# helper: run a command under /usr/bin/time -v, append its CPU (user+sys) to the
# running cpu accumulator file, and stream output to a per-stage console log.
CPU_ACC="logs/.cpu_acc"; : > "${CPU_ACC}"
timed() {  # <raw_out> <console_log> -- rest of args = command
  local raw="$1" clog="$2"; shift 2
  /usr/bin/time -v -o "${raw}" "$@" 2>&1 | tee "${clog}"
  awk -F': ' '/User time/{u=$2} /System time/{s=$2} END{printf "%.3f\n", u+s}' \
      "${raw}" >> "${CPU_ACC}"
}

# ---- STAGE 1: scm = base fit (node-0) + forward/backward selection ----------
echo "===== STAGE 1: scm base+search (${SCM_CFG}) ====="
# Stale-directory guard (resume after a killed run, e.g. TERM_RUNLIMIT): PsN's
# `scm` ABORTS if -directory=scm_dir already exists, and a partial scm_dir left
# by a walltime-killed attempt would make every resubmit of this cell fail
# instantly.  A cell is only "done" once logs/timing.json exists (written at the
# very end); if it does NOT, any scm_dir/refit here is an incomplete carcass ->
# wipe it so this run starts clean.  (SKIP_DONE at the array level already skips
# cells that DID finish, so this only ever fires on genuinely unfinished cells.)
if [ ! -f logs/timing.json ]; then
  rm -rf scm_dir refit
fi
timed logs/time_scm_raw.txt logs/scm_console.log \
  scm -config_file="${SCM_CFG}" -directory=scm_dir -nm_version="${NM_VERSION}" \
      ${RUN_ON_LSF_FLAG} -threads="${THREADS}" -clean=0
cp -f scm_dir/scmlog.txt logs/scmlog.txt 2>/dev/null || true

# ---- STAGE 2 (mandatory): final covariance refit of the winning model ------
# NOTE (PsN 5.5 on DVCHBS01): `scm` delegates "Writing final models" and
# `update_inits` to pharmpy, whose binary is MISSING in this build
# (".../bin/pharmpy: No such file or directory").  Consequence:
#   * scm_dir/final_models/final_backward.mod is written EMPTY/corrupt
#   * update_inits fails
# So we do NOT trust final_models/ and we do NOT use update_inits.  Instead we
# pick the REAL winning candidate model (a genuine NONMEM run written by PsN
# core) whose covariate-definition set matches the final selected set parsed
# from scmlog.txt, sanitize its header, append $COVARIANCE, and execute it.
# The refit is best-effort: if it cannot run, we WARN and continue (the parser
# still gets estimates + selection from the search; cov_done=FALSE).
echo "===== STAGE 2: final refit (execute winner + \$COVARIANCE) ====="
refit_ok=0

# final selected (param,covar) tags WITH state from scmlog.  PsN 5.5 writes each
# accepted step as "Parameter-covariate relation chosen in this forward step:
# CL-CRCL-4" (backward "... backward step: CL-SEX-1", state 1 == removed).  We
# keep the LAST forward state per (param,covar) minus backward-removed ones, so
# the signature is "PARAMCOVAR-STATE" (e.g. CLCRCL-5).  Matching on state as well
# as name is essential: a forward candidate can share the same covariate NAMES as
# the backward winner but hold one at a different state (exp/4 vs power/5).
kept_sig="$(awk '
  /relation chosen in this forward step:/  {n=split($NF,a,"-"); st[a[1] a[2]]=a[n]}
  /relation chosen in this backward step:/ {n=split($NF,a,"-");
                                            if(a[n]==1) delete st[a[1] a[2]];
                                            else        st[a[1] a[2]]=a[n]}
  END{ for(t in st) print t"-"st[t] }
' logs/scmlog.txt 2>/dev/null | sort | tr '\n' ' ' | sed 's/ *$//' || true)"
echo "final selected (tag-state): ${kept_sig}"

# tag+state signature of a candidate .mod (power=5, exp=4, cat=2).
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

# candidate models = real NONMEM runs under scm_dir (exclude PsN internals and
# the corrupt final_models); keep only files that start (after comments) with
# $PROB/$SIZES and whose tag-state signature exactly equals the selected set.
winner=""
while IFS= read -r m; do
  [ -f "${m}" ] || continue
  first="$(grep -vE '^\s*(;|$)' "${m}" | head -n1 || true)"
  case "${first}" in \$PROB*|\$SIZES*) : ;; *) continue ;; esac
  if [ "$(mod_sig "${m}")" = "${kept_sig}" ]; then winner="${m}"; break; fi
done < <(find scm_dir -name '*.mod' \
             -not -path '*/NM_run*' -not -path '*/temp_dir*' \
             -not -path '*/final_models/*' 2>/dev/null)

# NO base fallback: refitting the wrong (base or wrong-state) model would corrupt
# the covariate SE + condition number.  If no exact match, skip the refit.
if [ -z "${winner}" ]; then
  echo "WARN: no exact tag-state match for winner; skipping refit." >&2
fi

if [ -z "${winner}" ] || [ ! -f "${winner}" ]; then
  echo "WARN: no valid winning model found; skipping refit (cov_done will be FALSE)." >&2
else
  echo "winner model: ${winner}"
  # WIPE any stale refit artifacts first: PsN `execute -directory=refit_run
  # -clean=0` REUSES an existing refit_run/NM_run1 and copies its OLD results
  # back without re-running, silently yielding a base-model .lst/.ext (wrong OFV).
  rm -rf refit
  mkdir -p refit
  # Seed the refit from the winner's .ext final estimates (update_inits
  # replacement for the missing pharmpy) and force MAXEVAL=0: NONMEM then
  # computes the covariance AT the winner and the refit OFV is IDENTICAL to the
  # SCM winner OFV.  Executing the winner .mod as-is would re-minimise from
  # screening initials (covariate thetas ~0.001) and drift to a WRONG optimum,
  # giving SE / condition number for the wrong point in parameter space.
  winner_ext="${winner%.mod}.ext"
  seeded=0
  # ODE refit tol (structure=ode): export writes refit_tol.txt (TOL=/ATOL=); pass
  # it to seed_refit so the MAXEVAL=0 refit runs at the FINAL tol.  Absent for
  # advan4 (analytic) -> refit_tol_args stays empty and the winner .mod is used
  # as-is.
  refit_tol_args=()
  if [ -f "./refit_tol.txt" ]; then
    rtol_val="$(sed -n 's/^TOL=\([0-9]\+\).*/\1/p'  ./refit_tol.txt | head -1)"
    ratol_val="$(sed -n 's/^ATOL=\([0-9]\+\).*/\1/p' ./refit_tol.txt | head -1)"
    [ -n "${rtol_val}"  ] && refit_tol_args+=(--tol  "${rtol_val}")
    [ -n "${ratol_val}" ] && refit_tol_args+=(--atol "${ratol_val}")
  fi
  if [ -f "${winner_ext}" ]; then
    if run_seed --winner_mod "${winner}" --winner_ext "${winner_ext}" \
         --out refit/refit.mod --data "../data.csv" \
         --cov_file "${REFIT_COV_FILE}" --maxeval 0 "${refit_tol_args[@]}"; then
      seeded=1
    else
      echo "WARN: seed_refit failed; refit skipped (cov_done will be FALSE)." >&2
    fi
  else
    echo "WARN: winner .ext missing (${winner_ext}); refit skipped." >&2
  fi

  if [ "${seeded}" = "1" ]; then
    if ( cd refit && /usr/bin/time -v -o ../logs/time_refit_raw.txt \
           execute refit.mod -directory=refit_run -nm_version="${NM_VERSION}" \
                   ${RUN_ON_LSF_FLAG} -clean=0 ) 2>&1 | tee logs/refit_console.log
    then
      refit_ok=1
      awk -F': ' '/User time/{u=$2} /System time/{s=$2} END{printf "%.3f\n", u+s}' \
          logs/time_refit_raw.txt >> "${CPU_ACC}"
      cp -f refit/refit.lst logs/refit.lst 2>/dev/null || cp -f refit/refit_run/refit.lst logs/refit.lst 2>/dev/null || true
      cp -f refit/refit.ext logs/refit.ext 2>/dev/null || cp -f refit/refit_run/refit.ext logs/refit.ext 2>/dev/null || true
    else
      echo "WARN: covariance refit failed; continuing (cov_done will be FALSE)." >&2
    fi
  fi
fi

wall_end="$(date +%s.%N)"
end_iso="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# ---- totals ----------------------------------------------------------------
wall_sec="$(awk -v a="${wall_start}" -v b="${wall_end}" 'BEGIN{printf "%.3f", b-a}')"
cpu_sec="$(awk '{s+=$1} END{printf "%.3f", s}' "${CPU_ACC}")"
hog="$(awk -v c="${cpu_sec}" -v w="${wall_sec}" 'BEGIN{printf "%.3f",(w>0)?c/w:0}')"
scope="$([ "${refit_ok}" = "1" ] && echo scm+refit || echo scm)"

cat > logs/timing.json <<JSON
{
  "wall_sec": ${wall_sec},
  "cpu_sec": ${cpu_sec},
  "hog_factor": ${hog},
  "scope": "${scope}",
  "refit_ok": ${refit_ok},
  "threads": ${THREADS},
  "cores_per_fit": 1,
  "run_on_lsf": 0,
  "started": "${start_iso}",
  "ended": "${end_iso}",
  "nm_version": "${NM_VERSION}",
  "note": "Single-block total; scm(base+search)+refit; mirrors nlmixr2 t_base+t_scm+t_refit"
}
JSON
rm -f "${CPU_ACC}"

echo "TOTAL (${scope})  wall_sec=${wall_sec}  cpu_sec=${cpu_sec}  hog=${hog}  -> logs/timing.json"
