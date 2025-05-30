#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --tasks-per-node=8
#SBATCH -o TP1/output.txt
#SBATCH -e TP1/errors.txt
mpirun TP1/ejecutablepruebampi
