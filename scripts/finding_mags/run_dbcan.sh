#!/bin/bash

#######################################
# To predict CAZymes. More info can be found on our lab github (https://github.com/bcb-unl/run_dbcan)
# Input is protein sequences (here these are the ones predicted to contain coh/doc)
#######################################


input_dir="./uhgg_others/cow/protein_catalogue-100/putative_complex/faa_files"

base_output_dir="./uhgg_others/cow/protein_catalogue-100/putative_complex/dbcan_output"

mkdir -p "$base_output_dir"

# Loop through each .faa file in the input directory
for faa_file in "$input_dir"/*.faa; do
    # Extract the base name of the file without the directory and extension
    base_name=$(basename "$faa_file" .faa)
    
    # Define the output directory for each specific .faa file
    output_dir="${base_output_dir}/output_${base_name}"
    
    # Createe output directory for each faa
    mkdir -p "$output_dir"
    
    # Run dbCAN for the current .faa file, directing output to the dedicated directory
    run_dbcan "$faa_file" protein --out_dir "$output_dir"
    
    echo "Finished processing $faa_file"
done

echo "All files have been processed."
