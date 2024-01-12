#!/bin/sh
#SBATCH --job-name=raxml-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=6	#number of threads per task
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu

#raxml-ng-mpi --all --msa alignment_trimmed.phy --model TVM+I+G4 --prefix T3 --seed 2 --threads auto{11} --bs-trees 500 --bs-metric fbp,tbe

#raxml-ng-mpi --bootstrap --msa alignment_trimmed.phy --model TVM+I+G4 --prefix T7 --seed 7345 --threads 6 --bs-trees 500 --bs-metric fbp,tbe

raxml-ng-mpi --msa alignment_trimmed.phy --model TVM+I+G4 --prefix T9 --threads 6 --seed 224
