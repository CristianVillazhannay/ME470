#!/bin/bash
#SBATCH --job-name="poissonSolver"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=16  ## number of procs per node
#SBATCH -t 02:00:00
#SBATCH -o poissonSolver3D.out

module purge all
module load mpi 

mpirun -np 1 poissonSolver3D

mpirun -np 2 poissonSolver3D

mpirun -np 4 poissonSolver3D

mpirun -np 8 poissonSolver3D

mpirun -np 16 poissonSolver3D
