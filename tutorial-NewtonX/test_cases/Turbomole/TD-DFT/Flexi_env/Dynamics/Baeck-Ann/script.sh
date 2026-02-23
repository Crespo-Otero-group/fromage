#!/bin/bash -l
#
#$ -N TMDFTflBA
#$ -cwd # run in current working directlry
#$ -o ICOND_ADC # standard output file
#$ -l tmpfs=60G
#$ -pe smp 6 # number of cores and smp/parallel
#$ -l mem=8G # core memory (always 8GB for us)
#$ -l h_rt=11:59:0 # job time (240:0:0 or 1:0:0 for short queue)

export TURBODIR=/home/uccarcr/Programs/Turbomole_7p0/TURBOMOLE
export PARA_ARCH=SMP
export PATH=$PATH:$TURBODIR/scripts
export PATH=$PATH:$TURBODIR/bin/`sysname`
export PARNODES=$NSLOTS
export NX=/home/uccarcr/Programs/Newton-X/newtonx-cs/bin
module load perl/5.26.0
module load beta-modules
module load gcc-libs/10.2.0
module load compilers/intel/2022.2
export PATH=~/Programs/xtb-6.6.1/xtb/build:$PATH
module load python
source ~/Py-envs/fromage2p0/bin/activate

XTBPROCS=4

export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=6
export MKL_NUM_THREADS=6
export OMP_MAX_ACTIVE_LEVELS=1


ulimit -s unlimited

$NX/moldyn.pl > moldyn.log

