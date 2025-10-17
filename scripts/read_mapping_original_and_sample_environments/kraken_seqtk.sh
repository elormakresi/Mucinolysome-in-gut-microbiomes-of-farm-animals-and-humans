#!/bin/bash

## Quality control mostly for removal of contaminations (human reads)

# Define the path to your Kraken2 database
KRAKEN2_DB="/work/HCC/BCRF/app_specific/kraken/2.1"

# Define your input directory containing the FASTQ files
INPUT_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/bbduk_quality_processed"

# Define your output directory for storing non-human reads
OUTPUT_DIR="${INPUT_DIR}/human_reads_removed"
mkdir -p "$OUTPUT_DIR"

# Define the submission script file
submission_script="./read_files_putative_mucins/kraken_remove_human.sh"

# Initialize the submission script
echo "#!/bin/bash" > "$submission_script"

# Loop through all paired-end files and create individual SLURM scripts
for file1 in ${INPUT_DIR}/*_1.fastq; do
  # Get the base name without the _1.fastq suffix
  base=$(basename "$file1" _1.fastq)

  # Define file paths for paired-end reads
  file2="${INPUT_DIR}/${base}_2.fastq"

  # Define the SLURM script file name
  slurm_script="${OUTPUT_DIR}/run_kraken2_${base}.sh"

  # Create the SLURM script
  cat <<EOL > "$slurm_script"
#!/bin/bash
#SBATCH --time=40:00:00          # Run time in hh:mm:ss
#SBATCH --mem=200G                # Total memory required for the job
#SBATCH --cpus-per-task=1         # Request 2 cores
#SBATCH --job-name=kraken_${base}
#SBATCH --error=${OUTPUT_DIR}/kraken_${base}.%J.err
#SBATCH --output=${OUTPUT_DIR}/kraken_${base}.%J.out
#SBATCH --partition=yinlab,batch,guest

# Load the Kraken2 and seqtk modules
ml kraken2/2.1
ml seqtk/1.2

# Run Kraken2 to classify reads and filter out human reads
kraken2 --db "$KRAKEN2_DB" --paired "$file1" "$file2" \
  --output "${OUTPUT_DIR}/${base}_kraken2_output.txt" \
  --report "${OUTPUT_DIR}/${base}_kraken2_report.txt"

# Filter out only human reads and create a list of non-human read IDs (including unclassified reads)
# Modified to filter out host too
awk '\$3 != 9606 { print \$2 }' "${OUTPUT_DIR}/${base}_kraken2_output.txt" > "${OUTPUT_DIR}/${base}_non_human_read_ids.txt"
# awk '\$3 != 9598 && \$3 != 9443 { print \$2 }' "${OUTPUT_DIR}/${base}_kraken2_output.txt" > "${OUTPUT_DIR}/${base}_non_host_read_ids.txt"

# Extract non-human and host reads from the FastQ files
seqtk subseq "$file1" "${OUTPUT_DIR}/${base}_non_human_read_ids.txt" > "${OUTPUT_DIR}/${base}_1.fastq"
seqtk subseq "$file2" "${OUTPUT_DIR}/${base}_non_human_read_ids.txt" > "${OUTPUT_DIR}/${base}_2.fastq"

# Remove intermediate Kraken2 output files to save space
rm "${OUTPUT_DIR}/${base}_kraken2_output.txt" "${OUTPUT_DIR}/${base}_kraken2_report.txt" "${OUTPUT_DIR}/${base}_non_human_read_ids.txt"
EOL

  # Add the SLURM script to the submission script
  echo "sbatch $slurm_script" >> "$submission_script"
done

# Make the submission script executable
chmod +x "$submission_script"

echo "All SLURM scripts have been generated, and the submission script is ready."
echo "You can submit all jobs by running: ./$submission_script"
