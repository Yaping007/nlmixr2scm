#!/bin/bash
# ==============================================================================
# pilot_smoke.sh
# ------------------------------------------------------------------------------
# Pilot smoke test for the refit-true-model study.
#
#   3 cohorts (N40, N80, N300)
#   x 43 (scenario, boundary) combos per cohort
#   x 10 datasets per combo
#   = 1,290 fits
#
# At ~3 min average per fit and 50 parallel slots this finishes in ~90 min
# wall clock.  Use it to verify that:
#
#   1. All 129 job arrays queue without LSF syntax errors.
#   2. R / nlmixr2 / rxode2 load correctly on the compute nodes.
#   3. The master RDS files are readable from the compute nodes (paths).
#   4. Table-3 diagnostic fields populate for every boundary mode.
#   5. The PhysBnd rate for boundary in {none, wide} is >0 (expected).
#
# Once the pilot is green, re-submit with 250 datasets using
#   bash scripts/hpce/submit_all_arrays.sh 250 100
# ==============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 10 datasets per array, up to 50 array tasks running in parallel per array,
# all 3 cohorts (N40, N80, N300).
bash "$HERE/submit_all_arrays.sh" 10 50
