#!/bin/sh
#SBATCH --job-name=raxml-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=1	#number of threads per task
#SBATCH --time=35:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu


raxml-ng-mpi --all --msa alignment.phy --model TMV+I+G4 --prefix T3 --seed 2 --threads 1 --bs-trees 1000 --bs-metric fbp,tbe