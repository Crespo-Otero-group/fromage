#!/bin/sh
#
#$ -N Tr130_Eq0
#$ -cwd 
#$ -o output
#$ -j y
#$ -pe smp 4 
#$ -l h_vmem=8G 
#$ -l h_rt=240:0:0 

deactivate
source /data/home/btx721/Py-virtual-envs/fromage_devel/bin/activate
module load hdf5/1.10.2
module load intel/2018.3
module load lapack/3.9.0
export PATH=/data/SBCS-CrespoOtero/Fede/OpenMolcas-Intel/build:$PATH
export PATH=/data/SBCS-CrespoOtero/Fede/xtb/xtb-6.4.1/_build:$PATH # xtb for the QM' region
export MOLCAS_MEM=12000
export OMP_STACKSIZE=10G
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited

fro_run.py





