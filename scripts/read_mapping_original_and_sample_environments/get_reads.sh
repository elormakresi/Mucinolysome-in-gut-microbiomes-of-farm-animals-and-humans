#!/bin/bash

## to get read files from NCBI

# Load SRA toolkit module
ml SRAtoolkit/2.11

# Path to your list of SRA accession numbers
ACCESSIONS_FILE="./read_files_putative_mucins/testing_if_we_should_subsample/samples_to_test.txt"

# Base directory for downloaded files
BASE_DOWNLOAD_DIR="./read_files_putative_mucins/testing_if_we_should_subsample"

# Directory for SLURM scripts
SCRIPT_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/sample_download_scripts"
mkdir -p $SCRIPT_DIR

# Read each accession number and create a separate SLURM script
while read -r ACCESSION; do
    SCRIPT_FILE="$SCRIPT_DIR/download_$ACCESSION.sh"
    cat > $SCRIPT_FILE <<EOF
#!/bin/bash
#SBATCH --time=5:00:00          # Estimated time for download and processing
#SBATCH --mem-per-cpu=50G       # Memory per CPU
#SBATCH --job-name=download_$ACCESSION
#SBATCH --output=$ACCESSION.%J.out
#SBATCH --error=$ACCESSION.%J.err
#SBATCH --partition=yinlab,batch,guest

# Load SRA toolkit module
ml SRAtoolkit/2.11

# Set download directory
DOWNLOAD_DIR="$BASE_DOWNLOAD_DIR/$ACCESSION"
mkdir -p \$DOWNLOAD_DIR
cd \$DOWNLOAD_DIR

# Run fasterq-dump with --split-files for paired-end reads
fasterq-dump --split-files $ACCESSION

# Check if download was successful
if [ \$? -eq 0 ]; then
    echo "\$ACCESSION download completed successfully."
else
    echo "\$ACCESSION download failed." >> $BASE_DOWNLOAD_DIR/failed_downloads.txt
fi
EOF

    echo "SLURM script generated for $ACCESSION: $SCRIPT_FILE"
    # Optionally submit the job immediately
      sbatch $SCRIPT_FILE
done < "$ACCESSIONS_FILE"

echo "All scripts have been generated."
