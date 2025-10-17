#!/bin/bash

# Directory where the SLURM job scripts are stored
job_dir="/work/yinlab/jakresi/read_files_putative_mucins/testing_if_we_should_subsample/bbduk_slurm"

# Submit each job script
for job_script in "$job_dir"/*.sh; do
    sbatch "$job_script"
done

echo "All BBDuk jobs submitted."
