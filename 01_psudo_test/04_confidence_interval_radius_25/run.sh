#!/bin/bash
#SBATCH -p all
#SBATCH -J r25 # job name
#SBATCH -n 32  # num of total mpi processes
#SBATCH -c 1  # num of threads per mpi processes
#SBATCH -o run.log

# set GPU ID if needed
# export CUDA_VISIBLE_DEVICES="0"

# perform production with OpenMM
seq 100 | xargs -t -P10 -n1 ./run_each.sh

