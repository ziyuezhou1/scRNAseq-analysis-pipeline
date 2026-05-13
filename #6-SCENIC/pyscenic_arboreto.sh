#!/bin/bash

# SBATCH -p short    # partition name
#SBATCH --nodelist=cpu04
#SBATCH -c 30                   # Request N cores [>12 cores recommended for pyscenic2 steps]
#SBATCH -t 5-00:00               # hours:minutes time after which job will be killed [16h needed for 50k cells x 20k cells with 20 cores]
#SBATCH --mem 150G              # total amount of memory requested [>100G recommended]
#SBATCH --job-name henderson    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written


# Move to desired working directory
cd /share/home/hekwlab/ZZY/SCENIC

# Set up Python environment
# module load anaconda3
  ## This module is necessary to use *source* activate pyscenic2
# conda activate pyscenic

# Double-check what Python version is running
python3 --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH


# Create results directory
res_dir=outs_$(date +"%Y-%m-%d")
mkdir -p $res_dir
date

# Run pySCENIC get regulatory network from command line interface
# This circumvents the use of dask, which seems to cause trouble...
python3 ./arboreto_with_multiprocessing.py \
      /share/home/hekwlab/ZZY/SCENIC/so.renamed_res.0.2.loom \
        ./resources/allTFs_mm.txt \
        --num_workers 30 \
        -o $res_dir/pySCENIC_GRN_adjacencies.csv \
        --method grnboost2 \
        --sparse \
        --seed 787878

date
echo "regulatory networks inferred"