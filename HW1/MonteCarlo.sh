#!/bin/bash
#SBATCH --job-name="MonteCarlo"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=2  ## number of procs per node
#SBATCH -t 02:00:00
#SBATCH -o MonteCarlo.out

module purge all
module load mpi 
time mpirun -np 2 ./MonteCarlo 100000000

#!/bin/bash
#SBATCH --job-name="MonteCarlo"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=1  ## number of procs per node
#SBATCH -t 02:00:00
#SBATCH -o MonteCarlo.out

module purge all
module load mpi 
time mpirun -np 1 ./MonteCarlo 100000000

