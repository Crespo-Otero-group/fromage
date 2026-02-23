#!/bin/bash

#SBATCH --job-name=Orcadyn
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=0-0:59:00
#SBATCH --mem=40GB
#SBATCH --array=1-5

module load gcc/14.1.0-dn2g
module load openmpi/4.1.8-6env
export ORCAHOME=/user/home/Programs/Orca6.1.1/orca_6_1_1/
OPENMPI_PATH=/software/spack/linux-rocky8-broadwell/gcc-12.3.0/openmpi-4.1.8-6env/
export PATH=$OPENMPI_PATH/bin:$PATH
export LD_LIBRARY_PATH=$OPENMPI_PATH/lib:$PATH

export PATH=/user/home/Programs/dftbplus/install/arch/gcc14.1/bin:$PATH

source /user/home/Programs/anaconda3/etc/profile.d/conda.sh
conda activate fromage_devel

export NX=/user/home/Programs/newtonx-cs/bin
ulimit -s unlimited

TASK_DIR="TRAJ${SLURM_ARRAY_TASK_ID}"
cd $TASK_DIR/
$NX/moldyn.pl > moldyn.log

