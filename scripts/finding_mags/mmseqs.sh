#!/bin/bash

#######################################
#uses the human MAG  MGYG000005496's scaffoldins to find genomes. It was a filter instead of searching for coh/doc in  millions of MAGs
#queryDB has to be made  first before running this script. 
#######################################


# Navigate to the directory containing your .faa files
cd ./samodha_mags/s1
# Step 1: Create MMseqs2 database from .faa files
mmseqs createdb ./*.faa  targetDB

# Step 2: Create index for the database
mmseqs createindex targetDB tmp -s 7 --threads 32

# Step 3: Perform the search
mmseqs search ../../queryDB targetDB resultDB tmp --min-seq-id 0.50 --cov-mode 2 -c 0.85 -s 7 --threads 32

# Step 4: Convert alignment results to text format
mmseqs convertalis ../../queryDB targetDB resultDB mmseqs_results_samodha_s1.txt

echo "MMseqs2 workflow completed."

#Next part, I run on the command line
awk 'BEGIN{FS="\t"} $3 > 0.5 && $11 < 1e-5' mmseqs_results.txt > filtered_results.txt
