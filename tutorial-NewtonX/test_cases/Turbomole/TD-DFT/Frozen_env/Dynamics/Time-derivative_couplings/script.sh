#!/bin/bash -l
#
#$ -N TM-CIo_1
#$ -cwd # run in current working directlry
#$ -o NX_Gau # standard output file
#$ -l tmpfs=10G
#$ -pe smp 6 # number of cores and smp/parallel
#$ -l mem=8G # core memory (always 8GB for us)
#$ -l h_rt=47:59:0 # job time (240:0:0 or 1:0:0 for short queue)

export TURBODIR=/home/uccarcr/Programs/Turbomole_7p0/TURBOMOLE
export PARA_ARCH=SMP
export PATH=$PATH:$TURBODIR/scripts
export PATH=$PATH:$TURBODIR/bin/`sysname`
export PARNODES=$NSLOTS
export PATH=/home/uccarcr/Programs/xtb-6.5.1/build:$PATH
export NX=/home/uccarcr/Programs/Newton-X/newtonx-cs/bin
module load beta-modules
module load gcc-libs/10.2.0
module load perl/5.22.0
module load compilers/intel/2022.2
module load python
module load gaussian
source ~/Py-envs/fromage2p0/bin/activate

XTBPROCS=2

export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OMP_MAX_ACTIVE_LEVELS=

ulimit -s unlimited

$NX/moldyn.pl > moldyn.log

