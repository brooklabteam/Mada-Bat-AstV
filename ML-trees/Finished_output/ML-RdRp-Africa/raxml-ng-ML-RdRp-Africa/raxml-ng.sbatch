#!/bin/sh
#SBATCH --job-name=raxml-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=2	#number of threads per task
#SBATCH --time=30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu

#just to make sure it's loaded
#module load openmpi/3.0.0+gcc-10.1.0
#export PATH=/project2/cbrook/software/raxml-ng/bin:$PATH


#raxml-ng-mpi --all --msa alignment.phy --model GTR+I+G4 --prefix T3 --seed 2 --threads 2 --bs-trees 500 --bs-metric fbp,tbe

raxml-ng-mpi --bootstrap --msa alignment.phy --model GTR+I+G4 --prefix T7 --seed 222 --threads 2 --bs-trees 250 --bs-metric fbp,tbe
