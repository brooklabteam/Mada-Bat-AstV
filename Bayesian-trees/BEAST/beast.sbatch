#!/bin/sh
#SBATCH --job-name=beast
#SBATCH --account=pi-cbrook
#SBATCH --partition=broadwl
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8	#number of tasks
#SBATCH --cpus-per-task=1	#number of threads per task
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shorigan@uchicago.edu

module load beagle/trunk
module load beast/2.6.2

/software/beast-2.6.2-el7-x86_64/bin/beast -beagle_cpu -seed 777 astro_full_BEAST_alignment.xml 

#/software/beast-2.6.2-el7-x86_64/bin/beast -resume -threads 1 astro_full_BEAST_alignment.xml 
