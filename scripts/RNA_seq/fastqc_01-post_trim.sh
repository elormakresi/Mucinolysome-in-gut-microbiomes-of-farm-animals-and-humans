#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --mem=15gb
#SBATCH --job-name=fastqc_post_trim
#SBATCH --error=fastqc_post_trim.err
#SBATCH --output=fastqc_post_trim.out

## quality control
# Load necessary modules
ml fastqc/0.12
ml multiqc/py37/1.8

# Define the input directory for trimmed files
TRIMMED_DATA_DIR="./01.RawData/trim"

# Define a new output directory for post-trimming FastQC results
FASTQC_POST_TRIM_OUTPUT="./01.RawData/fastqc_output_post_trim"

# Create the output directory if it doesn't exist
mkdir -p "${FASTQC_POST_TRIM_OUTPUT}"

# Run FastQC on the trimmed paired-end files
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/G1_1_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/G1_2_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/G2_1_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/G2_2_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/G3_1_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/G3_2_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/M1_1_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/M1_2_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/M2_1_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/M2_2_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/M3_1_p.fq
fastqc -o "${FASTQC_POST_TRIM_OUTPUT}" "${TRIMMED_DATA_DIR}"/M3_2_p.fq


multiqc "${FASTQC_POST_TRIM_OUTPUT}"
