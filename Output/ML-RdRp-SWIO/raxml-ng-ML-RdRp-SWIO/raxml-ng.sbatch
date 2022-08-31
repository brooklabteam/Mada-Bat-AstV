#!/bin/sh
#SBATCH --job-name=raxml-ng-ML-RdRp-SWIO
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=1	#number of threads per task
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu


#raxml-ng-mpi --all --msa alignment.phy --model TVM+I+G4 --prefix T3 --seed 33 --threads 1 --bs-trees 500 --bs-metric fbp,tbe

#raxml-ng-mpi --msa alignment.phy --model TVM+I+G4 --prefix T5 --seed 2 --threads 1

raxml-ng-mpi --bootstrap --msa alignment.phy --model TVM+I+G4 --prefix T6 --seed 333 --threads 1 --bs-trees 500
