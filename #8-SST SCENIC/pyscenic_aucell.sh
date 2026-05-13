#!/bin/bash

#SBATCH --nodelist=cpu04
#SBATCH -c 12                   # Request N cores (must match future plan requested in R script)
#SBATCH -t 0-10:00               # hours:minutes runlimit after which job will be killed
#SBATCH --mem 100G               # total amount of memory requested
#SBATCH --job-name henderson    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written

# Double-check what Python version is running
python --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH


# Calculate cell module scores (from filtered loom, to extract acurate regulon incidence matrix later on)
pyscenic aucell \
    /share/home/hekwlab/ZZY/SST_SCENIC/so.renamed_res.0.2.loom \
    ./outs_2023-04-12/pySCENIC_CTX_regulons.csv \
    --output ./outs_2023-04-12/pySCENIC_filtered.loom \
    --num_workers 12

echo "cell scoring completed"

# Add dimensionality reduction (derived from pySCENIC AUCell matrix)
   ## I always skipped this step, as it can also be performed in R (with proper sample integration)
python add_visualization.py \
   --loom_input ./outs_2023-04-12/pySCENIC_filtered.loom \
   --loom_output ./outs_2023-04-12/pySCENIC_viz.loom \
   --num_workers 12

echo "dimensionality reduction added"