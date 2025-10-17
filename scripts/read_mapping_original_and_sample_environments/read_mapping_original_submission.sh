#!/bin/bash

# Directory where the SLURM scripts are stored
SLURM_SCRIPTS_DIR="/work/yinlab/jakresi/read_files_putative_mucins/ethiopian_MAG_samples/slurm_scripts"

# Loop through each SLURM script in the directory and submit it
for SLURM_SCRIPT in "$SLURM_SCRIPTS_DIR"/*.sh; do
    sbatch "$SLURM_SCRIPT"
done
