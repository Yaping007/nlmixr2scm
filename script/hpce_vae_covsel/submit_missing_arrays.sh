#!/bin/bash
# ============================================================================
# submit_missing_arrays.sh   (VAE covsel RESCUE / resubmit-gaps)
# ----------------------------------------------------------------------------
# Some array tasks exit without leaving a result (random node failures, OOM,
# walltime, silent VAE errors). This script scans the OUTPUT tree, works out
# which datasets in each (N, scenario, structure) cell are MISSING a result,
# and resubmits ONLY those dataset indices as a *discrete-index* LSF array
# (e.g.  vaecov_N300_scn16_linCmt_rescue[2,7,13-15]%40 ).
#
# A dataset id D in a cell is considered UNRUN when, under
#   <OUT_ROOT>/N<N>/scn<SS>_<structure>/covsel/
# EITHER
#   * res_ds<DDD>.rds is absent, OR
#   * res_ds<DDD>.rds is absent but res_ds<DDD>_ERROR.txt is present (failed).
# By default a present res_ds<DDD>.rds counts as DONE even if an
# _ERROR.txt also lingers; use REQUIRE_CLEAN=1 to also re-run any dataset
# that still has an _ERROR.txt sidecar (treat errored-but-wrote as unrun).
#
# It reuses the SAME task template (vae_covsel_array.lsf) and the SAME env-var
# contract as submit_one_array.sh, so rescued tasks are byte-for-byte the same
# computation as the original sweep — only the array index SET differs.
#
# ----------------------------------------------------------------------------
# Usage
#   bash script/hpce_vae_covsel/submit_missing_arrays.sh <n_ds> [maxpar] [ds_start]
#
#     n_ds     datasets-per-cell that the ORIGINAL sweep targeted (the range
#              checked is ds_start .. ds_start+n_ds-1)
#     maxpar   array throttle %N            (default 40)
#     ds_start first dataset index          (default 1)
#
#   Axis / tree overrides (env vars, space-separated where noted):
#     OUT_ROOT     output tree to scan   (default output/vae_covsel_pilot)
#     NS           "40 80 300"           (default 40 80 300)
#     SCENARIOS    "1 .. 16"             (default 1..16)
#     STRUCTURES   "linCmt ode"          (default linCmt ode)
#     REQUIRE_CLEAN=1   also re-run datasets that still carry an _ERROR.txt
#     DRYRUN=1     print the missing sets + bsub lines, submit NOTHING
#
# ----------------------------------------------------------------------------
# Examples
#   # Dry-run first: see exactly which datasets are missing in the 0723 run
#   OUT_ROOT=output/vae_covsel_full0723_est710 \
#   NS="40 80 300" SCENARIOS=16 STRUCTURES="linCmt ode" DRYRUN=1 \
#     bash script/hpce_vae_covsel/submit_missing_arrays.sh 245 40 1
#
#   # Then actually resubmit the gaps (drop DRYRUN):
#   OUT_ROOT=output/vae_covsel_full0723_est710 \
#   NS="40 80 300" SCENARIOS=16 STRUCTURES="linCmt ode" \
#     bash script/hpce_vae_covsel/submit_missing_arrays.sh 245 40 1
#
#   # Also re-run any dataset that left an _ERROR.txt (not just absent rds):
#   OUT_ROOT=output/vae_covsel_full0723_est710 REQUIRE_CLEAN=1 \
#     bash script/hpce_vae_covsel/submit_missing_arrays.sh 245 40 1
#
# Monitor:  bjobs -J 'vaecov_*_rescue'
# Kill:     bkill  -J 'vaecov_*_rescue'
# ============================================================================

set -euo pipefail

NDS=${1:?"Usage: $0 <n_ds> [maxpar] [ds_start]  (n_ds = datasets-per-cell the original sweep targeted)"}
MAXPAR=${2:-40}
DS_START=${3:-1}
DS_END=$((DS_START + NDS - 1))

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$HERE/../.." && pwd)}"

OUT_ROOT="${OUT_ROOT:-output/vae_covsel_pilot}"
REQUIRE_CLEAN="${REQUIRE_CLEAN:-0}"
DRYRUN="${DRYRUN:-0}"

[ -n "${NS:-}" ]         || NS="40 80 300"
[ -n "${SCENARIOS:-}" ]  || SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
[ -n "${STRUCTURES:-}" ] || STRUCTURES="linCmt ode"

# Auto-detect scripts subdir holding the driver (same logic as submit_one_array.sh)
SCRIPTS_DIR=""
for c in script scripts R; do
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

# ----------------------------------------------------------------------------
# Collapse a sorted list of integers into an LSF index spec: "2 3 4 7 9 10"
# -> "2-4,7,9-10".  LSF accepts commas + ranges inside the [ ] of -J.
# ----------------------------------------------------------------------------
collapse_ranges() {
  # reads space-separated sorted ints on stdin/args, echoes compact spec
  local nums=("$@") out="" start prev
  [ ${#nums[@]} -eq 0 ] && { echo ""; return; }
  start=${nums[0]}; prev=${nums[0]}
  local i n
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
echo "VAE covsel RESCUE scan"
echo "  OUT_ROOT      = $OUT_ROOT"
echo "  NS            = $NS"
echo "  SCENARIOS     = $SCENARIOS"
echo "  STRUCTURES    = $STRUCTURES"
echo "  ds range      = ${DS_START}..${DS_END}  (n_ds=$NDS)"
echo "  MAXPAR        = $MAXPAR"
echo "  REQUIRE_CLEAN = $REQUIRE_CLEAN   (1 = re-run datasets that left _ERROR.txt)"
echo "  DRYRUN        = $DRYRUN"
echo "=========================================================="

total_missing=0
total_cells=0
total_arrays=0

for N in $NS; do
  for SCN in $SCENARIOS; do
    SCN_PAD=$(printf '%02d' "$SCN")
    for ST in $STRUCTURES; do
      CELL_DIR="$REPO_ROOT/$OUT_ROOT/N${N}/scn${SCN_PAD}_${ST}/covsel"
      total_cells=$((total_cells + 1))

      # Build the list of missing dataset ids for this cell
      missing=()
      for ((D=DS_START; D<=DS_END; D++)); do
        DDD=$(printf '%03d' "$D")
        rds="$CELL_DIR/res_ds${DDD}.rds"
        err="$CELL_DIR/res_ds${DDD}_ERROR.txt"
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
      JOBNAME="vaecov_N${N}_scn${SCN_PAD}_${ST}_rescue"
      LOGDIR="$REPO_ROOT/logs/vaecov_N${N}"
      mkdir -p "$LOGDIR"

      echo "----------------------------------------------------------"
      echo "CELL  N${N} scn${SCN_PAD} ${ST}"
      echo "  dir     : $OUT_ROOT/N${N}/scn${SCN_PAD}_${ST}/covsel"
      echo "  missing : ${n_miss} dataset(s) -> [${IDX_SPEC}]"
      echo "  job     : ${JOBNAME}[${IDX_SPEC}]%${MAXPAR}"

      if [ "$DRYRUN" = "1" ]; then
        continue
      fi

      bsub \
        -J "${JOBNAME}[${IDX_SPEC}]%${MAXPAR}" \
        -o "${LOGDIR}/${JOBNAME}.%J.%I.out" \
        -e "${LOGDIR}/${JOBNAME}.%J.%I.err" \
        -env "all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, SAMPLE_N=${N}, SCN=${SCN}, STRUCTURE=${ST}, OUT_ROOT=${OUT_ROOT}" \
        < "$HERE/vae_covsel_array.lsf"
      total_arrays=$((total_arrays + 1))
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
  echo "Monitor : bjobs -J 'vaecov_*_rescue'"
  echo "Kill    : bkill  -J 'vaecov_*_rescue'"
fi
echo "=========================================================="
