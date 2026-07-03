#!/bin/bash
#BSUB -q short
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -M 2000
#BSUB -W 00:45

cd /home/liuya8j/nlmixr2scm/nlmixr2scm

export PATH=/CHBS/apps/EB/software/R/4.5.1-gomklc-2023b-0.2/bin:$PATH
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2

DS_ID=${DS_ID:-${LSB_JOBINDEX}}
SCN_ID=${SCN_ID:-16}
SAMPLE_N=${SAMPLE_N:-80}

echo "Running ds_id=$DS_ID scn=$SCN_ID N=$SAMPLE_N on $(hostname) at $(date)"

Rscript script/PerformanceEvaluation_sc_n_ds_HPCE.R $DS_ID $SCN_ID $SAMPLE_N

echo "Finished at $(date)"

