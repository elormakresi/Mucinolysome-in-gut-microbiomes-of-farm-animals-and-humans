#!/bin/bash

## may not be useful to you but our HPC allows for 1000 jobs to be submitted at a time so i used this to split the slurm scripts.

# Directory containing the .sh files
SOURCE_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/read_mapping_slurm_scripts_all"

# Base directory where the split folders will be created
TARGET_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/read_mapping_slurm_scripts_all"  

# Set the number of files per folder
FILES_PER_FOLDER=999

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Initialize counters
folder_counter=1
file_counter=0

# Create the first subdirectory
current_folder="${TARGET_DIR}/batch_${folder_counter}"
mkdir -p "$current_folder"

# Loop through the .sh files and move them in batches
for file in "${SOURCE_DIR}"/*.sh; do
    mv "$file" "$current_folder/"
    ((file_counter++))

    # When 900 files are moved, reset the counter and create a new folder
    if ((file_counter == FILES_PER_FOLDER)); then
        ((folder_counter++))
        current_folder="${TARGET_DIR}/batch_${folder_counter}"
        mkdir -p "$current_folder"
        file_counter=0
    fi
done

echo "Files have been split into batches of $FILES_PER_FOLDER."
