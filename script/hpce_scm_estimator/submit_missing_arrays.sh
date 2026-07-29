#!/bin/bash
# ============================================================================
# submit_missing_arrays.sh   (scm_bench estimator RESCUE / resubmit-gaps)
# ----------------------------------------------------------------------------
# Some array tasks exit without leaving a result (random node failures, OOM,
# walltime, silent fit errors). This script scans the OUTPUT tree, works out
# which datasets in each (N, scenario, structure, est, opt) cell are MISSING a
# result, and resubmits ONLY those dataset indices as a *discrete-index* LSF
# array (e.g.  bench_N300_scn16_linCmt_focei_bobyqa_rescue[2,7,13-15]%30 ).
#
# A dataset id D in a cell is UNRUN when, under
#   <OUT_ROOT>/N<N>/scn<SS>_<structure>/<est>_<opt>/
# EITHER
#   * res_ds<DDD>.rds is absent, OR
#   * res_ds<DDD>.rds is absent but <ds>_ERROR.txt is present (failed).
# A present res_ds<DDD>.rds counts as DONE even if an _ERROR.txt lingers;
# use REQUIRE_CLEAN=1 to also re-run any dataset that still carries an
# _ERROR.txt sidecar (treat errored-but-wrote as unrun).
#
# It reuses the SAME task template (bench_array.lsf) and the SAME env-var
# contract as submit_one_array.sh, so rescued tasks are byte-for-byte the same
# computation as the original sweep -- only the array index SET differs.
#
# ----------------------------------------------------------------------------
# Usage
#   bash script/hpce_scm_estimator/submit_missing_arrays.sh <n_ds> [maxpar] [ds_start]
#
#     n_ds     datasets-per-cell the ORIGINAL sweep targeted (range checked is
#              ds_start .. ds_start+n_ds-1)
#     maxpar   array throttle %N            (default 30)
#     ds_start first dataset index          (default 1)
#
#   Axis / tree overrides (env vars, space-separated where noted):
#     OUT_ROOT       output tree to scan    (default output/scm_bench)
#     NS             "40 80 300"            (default 40 80 300)
#     SCENARIOS      "1 .. 16"              (default 1..16)
#     STRUCTURES     "linCmt ode"           (default linCmt ode)
#     ESTIMATORS     "focei irlsfocei"      (default focei irlsfocei)
#     FOCEI_OPTS     "bobyqa"               (outer opts for focei row)
#     IRLSFOCEI_OPTS "bobyqa"               (outer opts for irlsfocei row)
#     FOCEIF_OPTS    "nlminb lbfgsb3c"      (parked foceif opts)
#     IRLS_OPTS      "lbfgsb3c"             (parked irlsfoceif opts)
#     REQUIRE_CLEAN=1   also re-run datasets that still carry an _ERROR.txt
#     DRYRUN=1          print the missing sets + bsub lines, submit NOTHING
#
#   The driver env knobs below are forwarded UNCHANGED (same defaults as
#   submit_one_array.sh) so rescued fits match the original sweep exactly:
#     WORKERS RX_THREADS SCREEN_SIGDIG SCREEN_ATOL SCREEN_RTOL WARM
#     REPACKAGE FORCE_RERUN R_LIBS_USER_OVERRIDE
#     QUEUE NCORES EXCLUSIVE      (LSF placement, as in submit_one_array.sh)
#
# ----------------------------------------------------------------------------
# Examples
#   # ifocei/bobyqa (the irlsfocei cell) ONLY, dry-run first:
#   OUT_ROOT=output/scm_bench \
#   ESTIMATORS=irlsfocei IRLSFOCEI_OPTS=bobyqa \
#   NS="40 80 300" SCENARIOS=16 STRUCTURES="linCmt ode" DRYRUN=1 \
#     bash script/hpce_scm_estimator/submit_missing_arrays.sh 250 30 1
#
#   # focei/bobyqa ONLY, actually resubmit the gaps:
#   OUT_ROOT=output/scm_bench \
#   ESTIMATORS=focei FOCEI_OPTS=bobyqa \
#   NS="40 80 300" SCENARIOS=16 STRUCTURES="linCmt ode" \
#     bash script/hpce_scm_estimator/submit_missing_arrays.sh 250 30 1
#
#   # both active cells (focei + irlsfocei, bobyqa), also re-run _ERROR.txt ones:
#   REQUIRE_CLEAN=1 \
#     bash script/hpce_scm_estimator/submit_missing_arrays.sh 250 30 1
#
# Monitor:  bjobs -J 'bench_*_rescue'
# Kill:     bkill  -J 'bench_*_rescue'
# ============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_ds> [maxpar] [ds_start]  (n_ds = datasets-per-cell the original sweep targeted)"}
MAXPAR=${2:-30}
DS_START=${3:-1}
DS_END=$((DS_START + NDS - 1))

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$HERE/../.." && pwd)}"

OUT_ROOT="${OUT_ROOT:-output/scm_bench}"
REQUIRE_CLEAN="${REQUIRE_CLEAN:-0}"
DRYRUN="${DRYRUN:-0}"

[ -n "${STRUCTURES:-}" ] || STRUCTURES="linCmt ode"
[ -n "${NS:-}" ]         || NS="40 80 300"
[ -n "${SCENARIOS:-}" ]  || SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
[ -n "${ESTIMATORS:-}" ] || ESTIMATORS="focei irlsfocei"
[ -n "${FOCEI_OPTS:-}" ]      || FOCEI_OPTS="bobyqa"
[ -n "${IRLSFOCEI_OPTS:-}" ]  || IRLSFOCEI_OPTS="bobyqa"
[ -n "${FOCEIF_OPTS:-}" ]     || FOCEIF_OPTS="nlminb lbfgsb3c"
[ -n "${IRLS_OPTS:-}" ]       || IRLS_OPTS="lbfgsb3c"

declare -A OPT_FOR
OPT_FOR[focei]="$FOCEI_OPTS"
OPT_FOR[foceif]="$FOCEIF_OPTS"
OPT_FOR[irlsfocei]="$IRLSFOCEI_OPTS"
OPT_FOR[irlsfoceif]="$IRLS_OPTS"

# Auto-detect scripts subdir holding the driver (same logic as submit_one_array.sh)
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

# Driver knobs forwarded unchanged (mirror submit_one_array.sh defaults) --------
WORKERS=${WORKERS:-3}
RX_THREADS=${RX_THREADS:-1}
SCREEN_SIGDIG=${SCREEN_SIGDIG:-NA}
SCREEN_ATOL=${SCREEN_ATOL:-NA}
SCREEN_RTOL=${SCREEN_RTOL:-NA}
WARM=${WARM:-calc}
REPACKAGE=${REPACKAGE:-0}
FORCE_RERUN=${FORCE_RERUN:-0}
R_LIBS_USER_OVERRIDE=${R_LIBS_USER_OVERRIDE:-$HOME/R/x86_64-pc-linux-gnu-library/4.5}

# ----------------------------------------------------------------------------
# Collapse a sorted list of ints into an LSF index spec:
#   "2 3 4 7 9 10" -> "2-4,7,9-10"  (LSF accepts commas + ranges inside [ ]).
# ----------------------------------------------------------------------------
collapse_ranges() {
  local nums=("$@") out="" start prev i n
  [ ${#nums[@]} -eq 0 ] && { echo ""; return; }
  start=${nums[0]}; prev=${nums[0]}
  for ((i=1; i<${#nums[@]}; i++)); do
    n=${nums[$i]}
    if [ "$n" -eq $((prev + 1)) ]; then
      prev=$n
    else
      if [ "$start" -eq "$prev" ]; then out+="${start},"; else out+="${start}-${prev},"; fi
      start=$n; prev=$n
    fi
  done
  if [ "$start" -eq "$prev" ]; then out+="${start}"; else out+="${start}-${prev}"; fi
  echo "$out"
}

echo "=========================================================="
echo "scm_bench estimator RESCUE scan"
echo "  OUT_ROOT      = $OUT_ROOT"
echo "  NS            = $NS"
echo "  SCENARIOS     = $SCENARIOS"
echo "  STRUCTURES    = $STRUCTURES"
echo "  ESTIMATORS    = $ESTIMATORS"
echo "  FOCEI_OPTS    = $FOCEI_OPTS"
echo "  IRLSFOCEI_OPTS= $IRLSFOCEI_OPTS"
echo "  ds range      = ${DS_START}..${DS_END}  (n_ds=$NDS)"
echo "  MAXPAR        = $MAXPAR"
echo "  REQUIRE_CLEAN = $REQUIRE_CLEAN   (1 = re-run datasets that left _ERROR.txt)"
echo "  DRYRUN        = $DRYRUN"
echo "=========================================================="

# Optional LSF placement knobs, exactly as submit_one_array.sh -----------------
EXTRA_BSUB=()
if [ -n "${QUEUE:-}" ]; then
  EXTRA_BSUB+=( -q "${QUEUE}" )
fi
if [ -n "${NCORES:-}" ]; then
  EXTRA_BSUB+=( -n "${NCORES}" -R "span[hosts=1]" )
fi
if [ "${EXCLUSIVE:-0}" = "1" ]; then
  EXTRA_BSUB+=( -x )
fi

total_missing=0
total_cells=0
total_arrays=0

for STRUCTURE in $STRUCTURES; do
  for N in $NS; do
    for SCN in $SCENARIOS; do
      SCN_PAD=$(printf '%02d' "$SCN")
      for EST in $ESTIMATORS; do
        OPTS="${OPT_FOR[$EST]:-}"
        if [ -z "$OPTS" ]; then
          echo "WARN: no valid outer_opt for est=$EST -- skipping" >&2
          continue
        fi
        for OPT in $OPTS; do
          CELL_DIR="$REPO_ROOT/$OUT_ROOT/N${N}/scn${SCN_PAD}_${STRUCTURE}/${EST}_${OPT}"
          total_cells=$((total_cells + 1))

          # Build list of missing dataset ids for this cell
          missing=()
          for ((D=DS_START; D<=DS_END; D++)); do
            DDD=$(printf '%03d' "$D")
            rds="$CELL_DIR/res_ds${DDD}.rds"
            err="$CELL_DIR/ds${DDD}_ERROR.txt"
            if [ ! -f "$rds" ]; then
              missing+=("$D")
            elif [ "$REQUIRE_CLEAN" = "1" ] && [ -f "$err" ]; then
              missing+=("$D")
            fi
          done

          n_miss=${#missing[@]}
          [ "$n_miss" -eq 0 ] && continue

          total_missing=$((total_missing + n_miss))
          IDX_SPEC=$(collapse_ranges "${missing[@]}")
          JOBNAME="bench_N${N}_scn${SCN_PAD}_${STRUCTURE}_${EST}_${OPT}_rescue"
          LOGDIR="$REPO_ROOT/logs/bench_N${N}"
          mkdir -p "$LOGDIR"

          echo "----------------------------------------------------------"
          echo "CELL  N${N} scn${SCN_PAD} ${STRUCTURE} ${EST}_${OPT}"
          echo "  dir     : $OUT_ROOT/N${N}/scn${SCN_PAD}_${STRUCTURE}/${EST}_${OPT}"
          echo "  missing : ${n_miss} dataset(s) -> [${IDX_SPEC}]"
          echo "  job     : ${JOBNAME}[${IDX_SPEC}]%${MAXPAR}"

          if [ "$DRYRUN" = "1" ]; then
            continue
          fi

          bsub \
            "${EXTRA_BSUB[@]}" \
            -J "${JOBNAME}[${IDX_SPEC}]%${MAXPAR}" \
            -o "${LOGDIR}/${JOBNAME}.%J.%I.out" \
            -e "${LOGDIR}/${JOBNAME}.%J.%I.err" \
            -env "all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, SAMPLE_N=${N}, SCN=${SCN}, EST=${EST}, OPT=${OPT}, STRUCTURE=${STRUCTURE}, WORKERS=${WORKERS}, RX_THREADS=${RX_THREADS}, OUT_ROOT=${OUT_ROOT}, SCREEN_SIGDIG=${SCREEN_SIGDIG}, SCREEN_ATOL=${SCREEN_ATOL}, SCREEN_RTOL=${SCREEN_RTOL}, WARM=${WARM}, REPACKAGE=${REPACKAGE}, FORCE_RERUN=${FORCE_RERUN}, R_LIBS_USER_OVERRIDE=${R_LIBS_USER_OVERRIDE}" \
            < "$HERE/bench_array.lsf"
          total_arrays=$((total_arrays + 1))
        done
      done
    done
  done
done

echo "=========================================================="
echo "Cells scanned      : $total_cells"
echo "Missing datasets   : $total_missing"
if [ "$DRYRUN" = "1" ]; then
  echo "DRYRUN=1 -> submitted NOTHING. Re-run without DRYRUN to resubmit."
else
  echo "Rescue arrays sent : $total_arrays"
  echo "Monitor : bjobs -J 'bench_*_rescue'"
  echo "Kill    : bkill  -J 'bench_*_rescue'"
fi
echo "=========================================================="
