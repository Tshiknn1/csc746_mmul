#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J csc746_f25_evan_caplinger_sum
#SBATCH --mail-user=ecaplinger@sfsu.edu
#SBATCH --mail-type=ALL
#SBATCH -t 00:05:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
srun -n 1 -c 256 --cpu_bind=cores build/benchmark-basic
srun -n 1 -c 256 --cpu_bind=cores build/benchmark-basic
srun -n 1 -c 256 --cpu_bind=cores build/benchmark-blas