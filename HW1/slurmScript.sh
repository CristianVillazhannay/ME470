#!/bin/bash
#SBATCH --job-name="lammps_lj_example"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=1  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 1 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log log1.log

#!/bin/bash
#SBATCH --job-name="lammps_lj_example"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=2  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 2 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log log2.log

#!/bin/bash
#SBATCH --job-name="lammps_lj_example"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=4  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 4 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log log4.log

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
mpirun -np 8 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log log8.log

#!/bin/bash
#SBATCH --job-name="lammps_lj_example"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 1                 ## number of nodes
#SBATCH --ntasks-per-node=16  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 16 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log log16.log

#!/bin/bash
#SBATCH --job-name="lammps_lj_example"
#SBATCH -A e31958
#SBATCH -p short             ## partition
#SBATCH -N 2                 ## number of nodes
#SBATCH --ntasks-per-node=16  ## number of procs per node
#SBATCH -t 02:00:00

## cd $PBS_O_WORKDIR
module purge all
module load lammps
mpirun -np 32 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log log32.log


