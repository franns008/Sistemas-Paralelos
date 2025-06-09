#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH --tasks-per-node=4

#SBATCH -o outputNonBlocking.txt
#SBATCH -e outputNonBlocking.txt
#SBATCH --overcommit

mpirun non-blocking 


