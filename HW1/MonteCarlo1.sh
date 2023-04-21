#!/bin/bash
#SBATCH --job-name="MonteCarlo"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=16  ## number of procs per node
#SBATCH -t 02:00:00
#SBATCH -o MonteCarlo.out

module purge all
module load mpi 

time mpirun -np 1 ./MonteCarlo 100000000

time mpirun -np 2 ./MonteCarlo 100000000

time mpirun -np 4 ./MonteCarlo 100000000

time mpirun -np 8 ./MonteCarlo 100000000

time mpirun -np 16 ./MonteCarlo 100000000
