#!/bin/sh
#SBATCH --job-name=raxml-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1 	#number of tasks
#SBATCH --cpus-per-task=1	#number of threads per task
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu

#just to make sure it's loaded
#module load openmpi/3.0.0+gcc-10.1.0
#export PATH=/project2/cbrook/software/raxml-ng/bin:$PATH


#raxml-ng-mpi --all --msa SWIO_RdRp_aligned_trimmed.phy --model TVM+I+G4 --prefix T3 --seed 2 --threads 1 --bs-trees 800 --bs-metric fbp,tbe
raxml-ng-mpi --bootstrap --msa SWIO_RdRp_aligned_trimmed.phy --model TVM+I+G4 --prefix T4 --seed 666 --threads 1 --bs-trees 400
