#!/bin/sh
#PBS -N Ni-test
#PBS -l select=1:ncpus=32:mem=96GB
#PBS -q normal
#PBS -j oe
#PBS -l walltime=50:00:00
#PBS -P Personal
 
cd ${PBS_O_WORKDIR}
 
## Capture Number of Cores
nprocs=`cat $PBS_NODEFILE|wc -l`
 
## Input File
#inputfile=
###Check mpi
which mpirun
echo -e "\n"
module purge
ulimit -s unlimited
source ~/VASP/vasp_env.sh
source /app/apps/oneapi/2022.2.0/setvars.sh --force
mpirun vasp_std
