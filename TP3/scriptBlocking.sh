#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH --tasks-per-node=4

#SBATCH -o outputBlocking.txt
#SBATCH -e outputBlocking.txt
#SBATCH --overcommit


mpirun blocking

