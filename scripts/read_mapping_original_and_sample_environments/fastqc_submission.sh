#!/bin/bash

# Directory where SLURM job scripts are stored
job_dir="/work/yinlab/jakresi/read_files_putative_mucins/testing_if_we_should_subsample/fastqc_scripts"

# Submit each job script
for job_script in "$job_dir"/*.sh; do
    sbatch "$job_script"
done

echo "All FastQC jobs submitted."
