#!/bin/bash
#BSUB -J "scm_full_scn16[10-250]%125"
#BSUB -q short
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -M 2000 
#BSUB -W 00:45
#BSUB -o logs/scm_full_scn16_%J_%I.out
#BSUB -e logs/scm_full_scn16_%J_%I.err

cd /home/liuya8j/nlmixr2scm/nlmixr2scm

export PATH=/CHBS/apps/EB/software/R/4.5.1-gomklc-2023b-0.2/bin:$PATH
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2

DS_ID=${LSB_JOBINDEX}
echo "Running ds_id = $DS_ID on host $(hostname) at $(date)"

Rscript script/PerformanceEvaluation_1DS_HPCE.R $DS_ID

echo "Finished at $(date)"