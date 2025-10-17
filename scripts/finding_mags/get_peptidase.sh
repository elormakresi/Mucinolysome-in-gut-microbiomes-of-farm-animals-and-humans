#!/bin/bash

#######################################
# search for peptidase in genomes
# input is protein sequences of coh/doc
#######################################


FAA_DIR="./putative_mucins_fna_kol109_duplicates_removed/faa_files_putative_mucins"

DB_PATH="./MEROPS.dmnd"

OUTPUT_DIR="./putative_mucins_fna_kol109_duplicates_removed/faa_files_putative_mucins/diamond_peptidase_no_coverage_sep_16"
mkdir -p "$OUTPUT_DIR"

# Loop over each .faa file in the directory
for faa_file in "$FAA_DIR"/*.faa
do
    # Extract filename without extension for output naming
    filename=$(basename "$faa_file" .faa)

    # Run Diamond blastp
    diamond blastp \
        --db "$DB_PATH" \
        --query "$faa_file" \
        --out "$OUTPUT_DIR/${filename}_results.txt" \
        --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid qlen slen length qstart qend sstart send qcovhsp scovhsp evalue bitscore --header \
        --evalue 1e-10 \
        --threads 24

    echo "Processed $faa_file"
done

echo "All Diamond searches completed."

