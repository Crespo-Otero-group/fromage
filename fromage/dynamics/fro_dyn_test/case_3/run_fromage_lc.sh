#!/bin/bash

export INPUT=molcas
export WORKDIR=/scratch/lijingbai2009/NBD_fromage/fro_dyn/case_3

if [ -d "/srv/tmp" ]
then
 export LOCAL_TMP=/srv/tmp
else
 export LOCAL_TMP=/tmp
fi

export MOLCAS=/work/lopez/Molcas
export MOLCAS_NPROCS=$SLURM_NTASKS
export MOLCAS_MEM=2250
export MOLCAS_PRINT=3
export MOLCAS_PROJECT=$INPUT
export MOLCAS_WORKDIR=$LOCAL_TMP/$USER/$SLURM_JOB_ID
export OMP_NUM_THREADS=1
export MOLCAS=/work/lopez/Molcas
export PATH=$MOLCAS/bin:$PATH

export PATH=/work/lopez/Python-3.7.4/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/lopez/Python-3.7.4/lib

export PATH=$PATH:/work/lopez/xtb-6.3.3/bin

cd $WORKDIR
fro_run.py
