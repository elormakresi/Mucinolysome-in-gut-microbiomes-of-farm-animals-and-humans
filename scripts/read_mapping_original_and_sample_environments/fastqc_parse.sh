#!/bin/bash

# Base directory where the FASTQ directories (subsampled) are located
BASE_DIR="/work/yinlab/jakresi/read_files_putative_mucins/testing_if_we_should_subsample/bbduk_quality_processed/human_reads_removed"

# Directory to store high-quality FASTQ files
HIGH_QUALITY_DIR="${BASE_DIR}/high_quality"

# Directory containing all the .out files from FastQC
OUT_FILES_DIR="/work/yinlab/jakresi/read_files_putative_mucins/testing_if_we_should_subsample/fastqc_results"

# Create the high-quality directory if it doesn't exist
mkdir -p "$HIGH_QUALITY_DIR"

# Loop through each .out file in the OUT_FILES_DIR
for OUT_FILE in "$OUT_FILES_DIR"/*.out; do
    echo "Processing $OUT_FILE..."

    # Variable to store the last identified FASTQ ID
    fastq_id=""

    # Parse each .out file line by line
    while read -r line; do
        # Check if the line contains the name of the paired FASTQ files to extract the FASTQ ID
        if [[ "$line" == *"inflating"* ]]; then
            # Extract the FASTQ ID from the line (assuming the filename format)
            fastq_id=$(echo "$line" | grep -oP '(ERR|SRR)\d+')
            echo "Identified FASTQ ID: $fastq_id"
        fi

        # Check if the line indicates that both files passed 'Per base sequence quality' and 'Per sequence quality scores'
        if [[ "$line" == *"Both paired files"* && "$line" == *"passed 'Per base sequence quality' and 'Per sequence quality scores'"* ]]; then
            if [ -n "$fastq_id" ]; then
                # Full path to the directory containing the FASTQ files
                full_fastq_dir="${BASE_DIR}/${fastq_id}"

                # Debugging: Print the full path and directory being checked
                echo "Checking directory: $full_fastq_dir"

                # Check if the directory exists before moving
                if [ -d "$full_fastq_dir" ]; then
                    echo "Directory $full_fastq_dir exists. Moving it."
                    # Move the directory to the high-quality directory
                    mv "$full_fastq_dir" "$HIGH_QUALITY_DIR/"
                    echo "Moved $full_fastq_dir to $HIGH_QUALITY_DIR"
                else
                    echo "Directory $full_fastq_dir does not exist. Skipping."
                fi
                # Clear the fastq_id variable after processing
                fastq_id=""
            fi
        fi
    done < "$OUT_FILE"
done

echo "All high-quality folders have been moved to $HIGH_QUALITY_DIR."
