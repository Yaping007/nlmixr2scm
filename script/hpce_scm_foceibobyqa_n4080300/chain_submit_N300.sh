#!/bin/bash
# Chain scenarios 1-16 with LSF dependencies for N=300.
# Each scenario starts only after the previous one ends.

PREV_JOB=""

for scn in $(seq 1 16); do
  scn_pad=$(printf "%02d" $scn)

  DEPENDENCY=""
  if [ -n "$PREV_JOB" ]; then
    DEPENDENCY="-w ended($PREV_JOB)"
  fi

  RESULT=$(SCN_ID=$scn SAMPLE_N=300 bsub $DEPENDENCY \
    -n 6 -M 8000 -W 04:00 -R "span[hosts=1]" \
    -J "scm_full_sc${scn_pad}_N300[1-250]%125" \
    -o logs/scm_full_sc${scn_pad}_N300_%J_%I.out \
    -e logs/scm_full_sc${scn_pad}_N300_%J_%I.err \
    < script/submit_full_sc_ds_N.sh)

  PREV_JOB=$(echo "$RESULT" | awk '{print $2}' | tr -d '<>')
  echo "Scenario $scn (N=300) submitted as job $PREV_JOB"
done

echo ""
echo "All 16 scenarios submitted for N=300 with dependencies."
echo "Only scenario 1 will run immediately."
echo "Scenarios 2-16 will PEND until predecessor finishes."