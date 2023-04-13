#!/bin/bash
#SBATCH --job-name="lammps_lj_example"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=8  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 8 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj
