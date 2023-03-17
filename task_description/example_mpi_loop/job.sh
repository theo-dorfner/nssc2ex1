#!/bin/bash

## job name

#SBATCH --job-name=mpi_loop

## logfiles (stdout/stderr) %x=job-name %j=job-id

#SBATCH --output=stdout-%x.%j.log
#SBATCH --error=stderr-%x.%j.log

## resource requests 

#SBATCH --partition=nssc    # partition for 360.242 and 360.242
#SBATCH --nodes=1           # request one node
#SBATCH --ntasks=40         # request eight processes on this node
#SBATCH --cpus-per-task=1   # request one cpu for each of these processes
#SBATCH --time=00:00:30     # set time limit to 30 seconds

## load modules and compilation (still on the login node)

module load pmi/pmix-x86_64     # [P]rocess [M]anagement [I]nterface (required by MPI-Implementation)
module load mpi/openmpi-x86_64  # MPI implementation (including compiler-wrappers mpicc/mpic++)

mpic++ -std=c++17 main.cpp -o main

## submitting jobs (on the allocated resources)

# job: print the hostname 

srun hostname

# job: run the mpi-enabled executable for different number of processes

for n in {1,2,5,10,20,30,40}
do
   srun --mpi=pmix --ntasks=${n} ./main
done
