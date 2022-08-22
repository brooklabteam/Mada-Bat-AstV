#!/bin/bash
#SBATCH --job-name=batAstro
#SBATCH --account=cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00

module load python

conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i alignment.fasta -t ml -p 8

scp -r shorigan@midway2.rcc.uchicago.edu:/home/shorigan/modeltest-ng/alignment.fasta.out /Users/sophiahorigan/Documents/GitHub/Mada-Bat-Astro/modeltest-ng/

conda deactivate
