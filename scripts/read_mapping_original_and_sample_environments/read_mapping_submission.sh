#!/bin/bash

SLURM_SCRIPTS_DIR="/work/yinlab/jakresi/read_files_putative_mucins/testing_if_we_should_subsample/read_mapping_scripts_scaf"
# Loop through each SLURM script and submit it
for SLURM_SCRIPT in "$SLURM_SCRIPTS_DIR"/*.sh; do
    sbatch "$SLURM_SCRIPT"
done
