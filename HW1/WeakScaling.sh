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
mpirun -np 1 lmp_Quest_gcc -var x 1 -var y 1 -var z 1 -in in.lj -log weaklog1.log

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
mpirun -np 2 lmp_Quest_gcc -var x 1 -var y 1 -var z 2 -in in.lj -log weaklog2.log

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
mpirun -np 4 lmp_Quest_gcc -var x 1 -var y 2 -var z 2 -in in.lj -log weaklog4.log

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
mpirun -np 8 lmp_Quest_gcc -var x 2 -var y 2 -var z 2 -in in.lj -log weaklog8.log

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
mpirun -np 16 lmp_Quest_gcc -var x 4 -var y 2 -var z 2 -in in.lj -log weaklog16.log

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
mpirun -np 32 lmp_Quest_gcc -var x 4 -var y 4 -var z 2 -in in.lj -log weaklog32.log


