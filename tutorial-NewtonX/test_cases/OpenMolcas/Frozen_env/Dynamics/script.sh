#!/bin/bash -l
#
#$ -N Molcas_3
#$ -cwd # run in current working directlry
#$ -o NX_Gau # standard output file
#$ -l tmpfs=10G
#$ -pe smp 6 # number of cores and smp/parallel
#$ -l mem=8G # core memory (always 8GB for us)
#$ -l h_rt=47:59:0 # job time (240:0:0 or 1:0:0 for short queue)

export PATH=/home/uccarcr/Programs/xtb-6.5.1_Intel2018/build:$PATH
export NX=/home/uccarcr/Programs/Newton-X/newtonx-cs/bin
export PATH=/home/uccarcr/Programs/Molcas/OpenMolcas_23.10/build:$PATH
module load default-modules
module load superlu/5.2.1/intel-2015-update2
module load arpack-ng/3.4.0/intel-2015-update2
module load armadillo/7.400.3/intel-2015-update2
module load hdf/5-1.10.2/intel-2018
module load python/3.8.6
source ~/Py-envs/fromage2p0/bin/activate

XTBPROCS=2
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export OMP_MAX_ACTIVE_LEVELS=1


ulimit -s unlimited

rm moldyn.log  -rf DEBUG/ TEMP/ INFO_RESTART/ RESULTS/
$NX/moldyn.pl > moldyn.log

