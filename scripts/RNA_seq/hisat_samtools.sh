#!/usr/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=128gb
#SBATCH --ntasks=8
#SBATCH --job-name=read_mapping
#SBATCH --partition=yinlab,batch,guest
#SBATCH --error=hisat_mapping.err
#SBATCH --output=hisat_mapping.out

module load hisat2/2.2
module load samtools/1.9

hisat2-build ./01.RawData/reference_based_rna_seq/kol109.fna ./01.RawData/reference_based_rna_seq/kol109_index

# --- Define Paths  ---
# Directory where HISAT2 index files are located and reference FASTA)
IDXPTH="./01.RawData/reference_based_rna_seq"
# Directory where your trimmed FASTQ files are located
READPTH="./01.RawData/trim"
# Directory to store your aligned BAM files and related outputs
BAM_DIR="./01.RawData/reference_based_rna_seq/output_bam_aligned"

# Name of your HISAT2 index prefix
HISAT2_INDEX_PREFIX="kol109_index"
# Name of your reference genome FASTA file (needed for samtools view -bT)
REFERENCE_FASTA="kol109.fna"

# Create the output directory if it doesn't exist
mkdir -p "$BAM_DIR"

echo "Starting HISAT2 alignment and BAM processing..."

# Loop through all your R1 (forward) trimmed FASTQ files
for R1_FILE in "${READPTH}"/*_1_p.fq; do
    # Extract the base name (e.g., G1, G2, M1)
    BASE_NAME=$(basename "$R1_FILE" _1_p.fq)

    # Construct the paths for R2 (reverse) file, SAM, BAM, and stats files
    R2_FILE="${READPTH}/${BASE_NAME}_2_p.fq"
    OUTPUT_SAM="${BAM_DIR}/${BASE_NAME}.sam"
    OUTPUT_BAM_SORTED="${BAM_DIR}/${BASE_NAME}_sorted.bam"
    OUTPUT_BAI="${BAM_DIR}/${BASE_NAME}_sorted.bai"
    OUTPUT_STATS="${BAM_DIR}/${BASE_NAME}.stats"
    OUTPUT_HISAT2_SUMMARY="${BAM_DIR}/${BASE_NAME}_hisat2_summary.txt"

    echo "--- Processing sample: ${BASE_NAME} ---"

    # Run HISAT2 alignment
    hisat2 -x "${IDXPTH}/${HISAT2_INDEX_PREFIX}" \
           -1 "$R1_FILE" \
           -2 "$R2_FILE" \
           -S "$OUTPUT_SAM" \
           --summary-file "$OUTPUT_HISAT2_SUMMARY" \
           -p 8

    # Check if HISAT2 command was successful
    if [ $? -ne 0 ]; then
        echo "HISAT2 alignment failed for ${BASE_NAME}. Exiting."
        exit 1
    fi

    echo "Converting SAM to BAM for ${BASE_NAME}..."
    # Convert SAM to BAM, using the reference FASTA for the header
    samtools view -bT "${IDXPTH}/${REFERENCE_FASTA}" "$OUTPUT_SAM" > "${BAM_DIR}/${BASE_NAME}.bam"

    # Check if SAM to BAM conversion was successful
    if [ $? -ne 0 ]; then
        echo "SAM to BAM conversion failed for ${BASE_NAME}. Exiting."
        exit 1
    fi

    echo "Sorting BAM for ${BASE_NAME}..."
    # Sort the BAM file by coordinate
    samtools sort "${BAM_DIR}/${BASE_NAME}.bam" -o "$OUTPUT_BAM_SORTED"

    # Check if sorting was successful
    if [ $? -ne 0 ]; then
        echo "BAM sorting failed for ${BASE_NAME}. Exiting."
        exit 1
    fi

    echo "Indexing sorted BAM for ${BASE_NAME}..."
    # Index the sorted BAM file
    samtools index "$OUTPUT_BAM_SORTED" "$OUTPUT_BAI"

    # Check if indexing was successful
    if [ $? -ne 0 ]; then
        echo "BAM indexing failed for ${BASE_NAME}. Exiting."
        exit 1
    fi

    echo "Generating alignment stats for ${BASE_NAME}..."
    # Generate alignment statistics per reference sequence
    samtools idxstats "$OUTPUT_BAM_SORTED" > "$OUTPUT_STATS"

    # Check if idxstats was successful
    if [ $? -ne 0 ]; then
        echo "samtools idxstats failed for ${BASE_NAME}. Exiting."
        exit 1
    fi

    echo "Removing intermediate unsorted BAM and SAM files for ${BASE_NAME}..."
    # Remove the intermediate unsorted BAM and SAM file to save space
    rm "${BAM_DIR}/${BASE_NAME}.bam" "$OUTPUT_SAM"

    echo "Finished processing sample: ${BASE_NAME}"
    echo "" # Add a blank line for readability between samples

done

echo "All HISAT2 alignment and BAM processing steps completed."
