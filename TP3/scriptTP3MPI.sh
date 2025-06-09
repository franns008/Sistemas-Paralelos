#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH --tasks-per-node=8
#SBATCH -o outputTP3MPIFran.txt
#SBATCH -e errorsTP3MPIFran.txt

mpirun TP3/resolucionMPIV2 $1
