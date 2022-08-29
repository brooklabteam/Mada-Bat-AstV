Code for running raxml-ng-mpi on the cluster
Great tutorial : https://github.com/amkozlov/raxml-ng/wiki/Tutorial

** DUPLICATE THE MASTER DIRECTORY AND MODIFY FOR EACH RUN **

1. Load Raxml

module load openmpi/3.0.0+gcc-10.1.0
export PATH=/project2/cbrook/software/raxml-ng/bin:$PATH



2. Ensure raxml works with your alignment

raxml-ng-mpi --check --msa alignment.phy --model GTR+I+G4 --prefix T1



3. Check for recommended threads

raxml-ng-mpi --parse --msa alignment.phy --model GTR+I+G4 --prefix T2



4. Update sbatch script with recommended threads & model



5. Run sbatch script

sbatch ./raxml-ng.sbatch



Check convergence

raxml-ng-mpi --bsconverge --bs-trees allbootstraps --prefix T8 --seed 2 --threads 1



Run additional bootstraps
raxml-ng --bootstrap --msa alignment.phy --model MODEL --prefix T5 --seed 333 --threads 1 --bs-trees #TREES --bs-metric fbp,tbe

cat T3.raxml.bootstraps T5.raxml.bootstraps > all bootstraps



After convergence, generate support values
raxml-ng-mpi --support --tree T3.raxml.bestTree --bs-trees allbootstraps --prefix T13 --threads 1

