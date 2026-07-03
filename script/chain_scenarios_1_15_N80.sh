#!/bin/bash
# Chain scenarios 1-15 with LSF dependencies.
# Each scenario starts only after the previous one ends.

PREV_JOB=""

for scn in $(seq 1 15); do
  scn_pad=$(printf "%02d" $scn)
  
  DEPENDENCY=""
  if [ -n "$PREV_JOB" ]; then
    DEPENDENCY="-w ended($PREV_JOB)"
  fi
  
  RESULT=$(SCN_ID=$scn SAMPLE_N=80 bsub $DEPENDENCY \
    -J "scm_full_sc${scn_pad}_N80[1-250]%250" \
    -o logs/scm_full_sc${scn_pad}_N80_%J_%I.out \
    -e logs/scm_full_sc${scn_pad}_N80_%J_%I.err \
    < script/submit_full_sc_ds_N.sh)
  
  PREV_JOB=$(echo "$RESULT" | awk '{print $2}' | tr -d '<>')
  echo "Scenario $scn submitted as job $PREV_JOB"
done

echo ""
echo "All scenarios 1-15 submitted with dependencies."
echo "Only scenario 1 will run immediately."
echo "Scenarios 2-15 will PEND until predecessor finishes."