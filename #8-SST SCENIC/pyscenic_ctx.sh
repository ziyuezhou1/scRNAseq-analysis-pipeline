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
cd /share/home/hekwlab/ZZY/SST_SCENIC

# Set up Python environment
# module load anaconda3
  ## This module is necessary to use *source* activate pyscenic2
# conda activate pyscenic

# Double-check what Python version is running
python3 --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH

# Clean up modules using feather ranking databases
pyscenic ctx outs_2023-10-12/pySCENIC_GRN_adjacencies.csv \
    databases/*.feather \
    --annotations_fname resources/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
    --expression_mtx_fname /share/home/hekwlab/ZZY/SST_SCENIC/so_inh.renamed_res0.2.loom \
    --output outs_2023-10-12/pySCENIC_CTX_regulons.csv \
    --mask_dropouts \
    --num_workers 30

date
echo "regulon modules defined"
sleep 5