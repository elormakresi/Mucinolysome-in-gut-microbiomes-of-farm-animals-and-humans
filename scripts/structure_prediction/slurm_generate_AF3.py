#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 14 17:52:26 2025

@author: siva
"""

import os
from glob import glob

# Template SLURM script
template = """#!/bin/bash
#SBATCH --job-name=AlphaFold_{name}_GPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --partition=guest_gpu
#SBATCH --gres=gpu
#SBATCH --error=AlphaFold_{name}_GPU.%J.err
#SBATCH --output=AlphaFold_{name}_GPU.%J.out

module purge
module load alphafold3/3.0

apptainer exec docker://unlhcc/alphafold3:3.0.1 run_alphafold.sh \\
    --model_dir=$AF3_MODELS \\
    --db_dir=$AF3_DBS \\
    --json_path={name}.json \\
    --output_dir=/AF3_combo/{name}
"""

# Get all JSON files in current directory
json_files = glob("*.json")

for json_file in json_files:
    name = os.path.splitext(os.path.basename(json_file))[0]
    slurm_filename = f"run_{name}.slurm"

    with open(slurm_filename, "w") as f:
        f.write(template.format(name=name, json_file=os.path.abspath(json_file)))

    print(f"Generated: {slurm_filename}")
