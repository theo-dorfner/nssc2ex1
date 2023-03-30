#!/bin/bash

#SBATCH --job-name=jacobiMPI_1D

#SBATCH --output=stdout-%x%j.txt

## resource requests

#SBATCH --partition=nssc
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

## load modules

module load pmi/pmix-x86_64
module load mpi/openmpi-x86_64

## start compilation
mpic++ -std=c++17 -O3 -Wall -pedantic -march=native -ffast-math ../src/jacobiMPI.cpp -o jacobiMPI

## now this is the proper job
for n in {1,2,3,5,7,10,12,15,17,20,22,25,27,30,32,35,37,40}
do
    srun --mpi=pmix --ntasks=${n} ./jacobiMPI 1D 250 30
done