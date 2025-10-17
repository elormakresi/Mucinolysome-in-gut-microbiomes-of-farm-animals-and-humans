#!/bin/bash

#######################################
# to identify cohesins in protein sequences using hmmer
# input is cohesin.hmm and protein sequence file
#######################################


# Directory containing .faa files
DIR="./faa_files_putative_mucins/mucosal_ref_genomes"

# Output directory for .tsv files
OUTPUT_DIR="./faa_files_putative_mucins/mucosal_ref_genomes/coh_doc_results"

# Ensure the output directory exists
mkdir -p $OUTPUT_DIR

# Loop over each .faa file in the directory
for faa_file in $DIR/*.faa; do
    # Extract the base name for the output file
    base_name=$(basename $faa_file .faa)

    # Run hmmsearch for each file
    hmmsearch -o /dev/null --noali --cpu 32 -E 0.0001 --domE 0.0001 --domtblout "${OUTPUT_DIR}/${base_name}_cohesin_domtbl.tsv" cohesin.hmm $faa_file
done
