#!/bin/bash
# ==============================================================================
# submit_one_array.sh
# ------------------------------------------------------------------------------
# Submit ONE LSF job array for the refit-true-model study covering
#   cohort x scenario x boundary  x  {dataset 1 .. dataset N}
#
# The script auto-detects its own directory ("$HERE") and reads
# refit_array.lsf from the same directory, so it works regardless of whether
# you keep the scripts in scripts/hpce/, script/hpce_refit/, or anywhere else.
#
# Repo root defaults to $HERE/../.. (i.e. scripts live 2 levels below the
# repo root).  Override by exporting REPO_ROOT before invoking, e.g.
#   REPO_ROOT=/home/liuya8j/nlmixr2scm/nlmixr2scm bash <this>/submit_one_array.sh ...
#
# Usage:
#   bash <hpce_dir>/submit_one_array.sh <cohort> <scenario> <boundary> <n_datasets> [max_parallel]
#
# Example (pilot: 10 datasets, up to 50 tasks in parallel):
#   bash script/hpce_refit/submit_one_array.sh N300 16 none 10 50
#
# Example (full: 250 datasets):
#   bash script/hpce_refit/submit_one_array.sh N300 16 narrow 250 100
# ==============================================================================

set -euo pipefail

COHORT=${1:?"Usage: $0 <cohort:N40|N80|N300> <scenario:2..16> <boundary:none|wide|narrow|tight> <n_datasets> [max_parallel]"}
SCN=${2:?"scenario (2..16) required"}
BND=${3:?"boundary (none|wide|narrow|tight) required"}
NDS=${4:?"n_datasets required (e.g. 10 for pilot, 250 for full)"}
MAXPAR=${5:-50}

# ---- Resolve HPCE-script dir + repo root ---------------------------------
# HERE        = directory containing THIS script (always correct, any layout)
# REPO_ROOT   = workspace root; default = $HERE/../.., override via env var
# SCRIPTS_DIR = subdir of REPO_ROOT holding refit_one_dataset.R.  Auto-detected
#               by trying common names ("scripts", "script"), then falling
#               back to a repo-wide `find`.  Override via env var if needed:
#                   SCRIPTS_DIR=my_r_dir bash <this>/submit_one_array.sh ...
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$HERE/../.." && pwd)}"

# 1. Explicit override wins.
# 2. Otherwise try common names in the workspace layout.
# 3. Otherwise fall back to a repo-wide search (max depth 4).
if [ -n "${SCRIPTS_DIR:-}" ]; then
    if [ ! -f "$REPO_ROOT/$SCRIPTS_DIR/refit_one_dataset.R" ]; then
        echo "ERROR: SCRIPTS_DIR='$SCRIPTS_DIR' is set but" >&2
        echo "       $REPO_ROOT/$SCRIPTS_DIR/refit_one_dataset.R does not exist." >&2
        exit 1
    fi
else
    SCRIPTS_DIR=""
    for candidate in scripts script R inst/scripts; do
        if [ -f "$REPO_ROOT/$candidate/refit_one_dataset.R" ]; then
            SCRIPTS_DIR="$candidate"
            break
        fi
    done
    # Deep search fallback (relative path from REPO_ROOT, up to 4 levels deep)
    if [ -z "$SCRIPTS_DIR" ]; then
        found="$(find "$REPO_ROOT" -maxdepth 4 -name refit_one_dataset.R \
                     -not -path '*/.*' 2>/dev/null | head -n 1)"
        if [ -n "$found" ]; then
            SCRIPTS_DIR="$(dirname "${found#$REPO_ROOT/}")"
        fi
    fi
fi

# ---- Sanity checks --------------------------------------------------------
if [ ! -f "$HERE/refit_array.lsf" ]; then
    echo "ERROR: refit_array.lsf not found at $HERE/refit_array.lsf" >&2
    exit 1
fi
if [ -z "$SCRIPTS_DIR" ]; then
    echo "ERROR: refit_one_dataset.R not found anywhere under REPO_ROOT=$REPO_ROOT" >&2
    echo "" >&2
    echo "  Tried standard locations:" >&2
    for candidate in scripts script R inst/scripts; do
        echo "    $REPO_ROOT/$candidate/refit_one_dataset.R" >&2
    done
    echo "" >&2
    echo "  Deep search (find -maxdepth 4) also returned nothing." >&2
    echo "" >&2
    echo "  R files found under REPO_ROOT:" >&2
    find "$REPO_ROOT" -maxdepth 4 -name '*.R' -not -path '*/.*' 2>/dev/null \
         | head -n 20 | sed 's/^/    /' >&2
    echo "" >&2
    echo "  Fixes:" >&2
    echo "    1. Copy refit_one_dataset.R, refit_helpers.R, true_model_factory.R" >&2
    echo "       to a directory under REPO_ROOT (e.g. \$REPO_ROOT/script/)." >&2
    echo "    2. Or set SCRIPTS_DIR to the correct subdir before invoking:" >&2
    echo "         SCRIPTS_DIR=my_dir bash $0 ..." >&2
    exit 1
fi

# ---- Log directory (per cohort so subdirs stay tidy) -----------------------
LOGDIR="$REPO_ROOT/logs/refit_${COHORT}"
mkdir -p "$LOGDIR"

# ---- LSF job name (used in -J/-o/-e/tokill patterns) -----------------------
SCN_PAD=$(printf '%02d' "$SCN")
JOBNAME="refit_${COHORT}_scn${SCN_PAD}_${BND}"

# ---- Expand path templates ------------------------------------------------
# MASTER_RDS_TEMPLATE lets you use one env var across all 129 arrays.  Tokens:
#   {COHORT} -> "N40" / "N80" / "N300"
#   {SCN}    -> 2-digit zero-padded scenario, e.g. "04", "16"
# Example for the HPCE per-scenario layout:
#   export MASTER_RDS_TEMPLATE='Inputdataset/sim_obs_{COHORT}/sim_obs_scenario_{SCN}.rds'
#
# ---- HPCE auto-detection ---------------------------------------------------
# If neither MASTER_RDS nor MASTER_RDS_TEMPLATE was set by the caller AND the
# HPCE-conventional Inputdataset/ directory exists at REPO_ROOT, assume the
# HPCE per-scenario layout and set both MASTER_RDS_TEMPLATE and OUT_ROOT to
# the values that made the earlier "narrow" runs work.  This prevents the
# R driver from silently falling back to its Windows workspace defaults
# ("simulated_virtual_dataset_*/sim_obs_all_scenarios.rds" +
#  "outputs/refit_true_*/") on HPCE, which produce
#   ERROR: file.exists(master_rds) is not TRUE
# and outputs landing in the wrong directory tree.
if [ -z "${MASTER_RDS:-}" ] && [ -z "${MASTER_RDS_TEMPLATE:-}" ] \
        && [ -d "$REPO_ROOT/Inputdataset" ]; then
    MASTER_RDS_TEMPLATE='Inputdataset/sim_obs_{COHORT}/sim_obs_scenario_{SCN}.rds'
    : "${OUT_ROOT:=output}"
    echo "[hpce-auto] Inputdataset/ found -> using HPCE defaults:"
    echo "            MASTER_RDS_TEMPLATE='$MASTER_RDS_TEMPLATE'"
    echo "            OUT_ROOT='$OUT_ROOT'"
fi

if [ -z "${MASTER_RDS:-}" ] && [ -n "${MASTER_RDS_TEMPLATE:-}" ]; then
    MASTER_RDS="${MASTER_RDS_TEMPLATE//\{COHORT\}/$COHORT}"
    MASTER_RDS="${MASTER_RDS//\{SCN\}/$SCN_PAD}"
fi

# OUT_ROOT lets you redirect the output tree.  Full path becomes
#   $OUT_ROOT/refit_true_<COHORT>/scn<NN>_<bnd>/
# Default (unset): R driver builds outputs/refit_true_<COHORT>/scn<NN>_<bnd>/
if [ -z "${OUT_DIR:-}" ] && [ -n "${OUT_ROOT:-}" ]; then
    OUT_DIR="${OUT_ROOT}/refit_true_${COHORT}/scn${SCN_PAD}_${BND}"
fi

# ---- Sanity check the master RDS if we can resolve it ---------------------
if [ -n "${MASTER_RDS:-}" ]; then
    # If path is not absolute, resolve relative to REPO_ROOT
    if [ "${MASTER_RDS#/}" = "$MASTER_RDS" ]; then
        MASTER_RDS_ABS="$REPO_ROOT/$MASTER_RDS"
    else
        MASTER_RDS_ABS="$MASTER_RDS"
    fi
    if [ ! -f "$MASTER_RDS_ABS" ]; then
        echo "ERROR: master RDS not found at: $MASTER_RDS_ABS" >&2
        echo "       (derived from MASTER_RDS='$MASTER_RDS')" >&2
        echo "" >&2
        echo "  Fixes:" >&2
        echo "    1. Verify the file exists at that path." >&2
        echo "    2. Or set MASTER_RDS_TEMPLATE explicitly, e.g.:" >&2
        echo "         export MASTER_RDS_TEMPLATE='Inputdataset/sim_obs_{COHORT}/sim_obs_scenario_{SCN}.rds'" >&2
        exit 1
    fi
fi

echo "Submitting array: ${JOBNAME}[1-${NDS}]%${MAXPAR}"
echo "  HERE        = $HERE"
echo "  REPO_ROOT   = $REPO_ROOT"
echo "  SCRIPTS_DIR = $SCRIPTS_DIR    (contains refit_one_dataset.R)"
echo "  MASTER_RDS  = ${MASTER_RDS:-<cohort default>}"
echo "  OUT_DIR     = ${OUT_DIR:-<cohort default>}"
echo "  logs     -> $LOGDIR/${JOBNAME}.<jobid>.<idx>.{out,err}"

# Build the bsub -env string.  Optional overrides only appear if the caller
# exported them; otherwise the LSF template falls through to cohort defaults.
ENV_LIST="all, REPO_ROOT=${REPO_ROOT}, SCRIPTS_DIR=${SCRIPTS_DIR}, COHORT=${COHORT}, SCN=${SCN}, BND=${BND}"
if [ -n "${MASTER_RDS:-}" ]; then
    ENV_LIST="${ENV_LIST}, MASTER_RDS=${MASTER_RDS}"
fi
if [ -n "${OUT_DIR:-}" ]; then
    ENV_LIST="${ENV_LIST}, OUT_DIR=${OUT_DIR}"
fi

bsub \
    -J "${JOBNAME}[1-${NDS}]%${MAXPAR}" \
    -o "${LOGDIR}/${JOBNAME}.%J.%I.out" \
    -e "${LOGDIR}/${JOBNAME}.%J.%I.err" \
    -env "$ENV_LIST" \
    < "$HERE/refit_array.lsf"
