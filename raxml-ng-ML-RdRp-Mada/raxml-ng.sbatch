#!/bin/sh
#SBATCH --job-name=raxml-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=1	#number of threads per task
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu

#ML and bootstraps
#raxml-ng-mpi --all --msa alignment.phy --model TPM2uf+I+G4 --prefix T3 --seed 2 --threads 1 --bs-trees 1000 --bs-metric fbp,tbe

#bootstraps only
#use different seed to ensure independence
raxml-ng-mpi --bootstrap --msa alignment.phy --model TPM2uf+I+G4 --prefix T12 --seed 88 --threads 1 --bs-trees 100 --bs-metric fbp,tbe
