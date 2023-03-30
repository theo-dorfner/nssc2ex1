#!/bin/bash

#SBATCH --job-name=jacobiMPI_1D

#SBATCH --output=stdout-%x%j.txt

## resource requests

#SBATCH --partition=nssc
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00 #set time to 1h30min of intial tests

## load modules

module load pmi/pmix-x86_64
module load mpi/openmpi-x86_64

## start compilation (?should we compile before that? i.e. manually to catch errors)
mpic++ -std=c++17 -O3 -Wall -pedantic -march=native -ffast-math ../src/jacobiMPI.cpp -o jacobiMPI


## submitting jobs

## I don't know why they print the hostname in the original job submission, but yeah

##srun hostname

## now this is the proper job
## ./jacobiMPI <resolution> <iterations>

for n in {1,2,3,5,7,10,12,15,17,20,22,25,27,30,32,35,37,40}
do
    srun --mpi=pmix --ntasks=${n} ./jacobiMPI 1D 1000 30
done