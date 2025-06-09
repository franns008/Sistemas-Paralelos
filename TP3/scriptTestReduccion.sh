#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH --tasks-per-node=8

#SBATCH -o outputTestReduccion.txt
#SBATCH -e outputTestReduccion.txt
#SBATCH --overcommit


mpirun testReduccion $1

