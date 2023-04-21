#!/bin/bash
#SBATCH --job-name="weakproc32"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 2                 ## number of nodes
#SBATCH --ntasks-per-node=16  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 32 lmp_Quest_gcc -var x 4 -var y 4 -var z 2 -in in.lj -log weaklog32.log




