#!/bin/sh
#SBATCH --job-name=raxml-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=6	#number of threads per task
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu

raxml-ng-mpi --all --msa alignment.phy --model TVM+I+G4 --prefix T3 --seed 2 --threads 6 --bs-trees 1000 --bs-metric fbp,tbe
