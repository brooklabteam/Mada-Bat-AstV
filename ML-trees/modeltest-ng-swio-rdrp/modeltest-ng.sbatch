#!/bin/bash
#SBATCH --job-name=modeltest-ng
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00

module load python

conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i SWIO_RdRp_aligned_trimmed.fasta -d nt -t ml -p 8 

conda deactivate
