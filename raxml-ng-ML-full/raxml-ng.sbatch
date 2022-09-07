#!/bin/sh
#SBATCH --job-name=raxml-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=1	#number of threads per task
#SBATCH --time=05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu

#just to make sure it's loaded
module load openmpi/3.0.0+gcc-10.1.0
export PATH=/project2/cbrook/software/raxml-ng/bin:$PATH


raxml-ng-mpi --all --msa alignment.phy --model MODEL --prefix T3 --seed 2 --threads 1 --bs-trees 1000 --bs-metric fbp,tbe