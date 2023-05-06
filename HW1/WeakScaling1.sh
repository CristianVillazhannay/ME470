#!/bin/bash
#SBATCH --job-name="weakproc16"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=16  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 1 lmp_Quest_gcc -var x 1 -var y 1 -var z 1 -in in.lj -log weaklog1.log

mpirun -np 2 lmp_Quest_gcc -var x 1 -var y 1 -var z 2 -in in.lj -log weaklog2.log

mpirun -np 4 lmp_Quest_gcc -var x 1 -var y 2 -var z 2 -in in.lj -log weaklog4.log

mpirun -np 8 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log weaklog8.log

mpirun -np 16 lmp_Quest_gcc -var x 4 -var y 2 -var z 2 -in in.lj -log weaklog16.log


