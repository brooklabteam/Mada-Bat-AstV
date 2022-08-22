#!/bin/bash
#SBATCH --job-name=raxml-ng
#SBATCH --account=cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu


raxml-ng-mpi --all --msa alignment.phy --model TPM2uf+I+G4 --prefix T3 --seed 2 --threads 1 --bs-metric fbp,tbe