#!/bin/bash
#SBATCH --nodelist=cpu06
#SBATCH -n 80
Rscript generate_exprmat.R
date -R
# mpirun hostname
sleep 10
echo Completed!