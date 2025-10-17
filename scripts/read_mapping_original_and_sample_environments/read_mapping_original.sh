#!/bin/bash

# Directories
BASE_DIR="./read_files_putative_mucins/ethiopian_MAG_samples"
SLURM_SCRIPTS_DIR="$BASE_DIR/slurm_scripts"  

# Create the SLURM scripts directory if it does not exist
mkdir -p "$SLURM_SCRIPTS_DIR"

# Loop through each sample directory in the base directory
for SAMPLE_DIR in "$BASE_DIR"/*; do
    SAMPLE_ID=$(basename "$SAMPLE_DIR")
    GENOME_PATH=$(find "$SAMPLE_DIR" -name '*.fna')
    GENOME=$(basename "$GENOME_PATH" .fna)

    # Create SLURM script for this mapping
    SLURM_SCRIPT="${SLURM_SCRIPTS_DIR}/${SAMPLE_ID}_mapping.sh"
    cat <<EOF > "$SLURM_SCRIPT"
#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name=${SAMPLE_ID}_mapping
#SBATCH --partition=yinlab,batch,guest

ml bowtie/2.3
ml instrain/1.3
ml coverm/0.6
ml samtools/1.19

# Build index
bowtie2-build "$GENOME_PATH" "${SAMPLE_DIR}/${GENOME}_index"

# Check if paired-end or single-end
if [ -f "${SAMPLE_DIR}/${SAMPLE_ID}_2.fastq" ]; then
    # Paired-end mapping
    bowtie2 -x "${SAMPLE_DIR}/${GENOME}_index" \\
            -1 "${SAMPLE_DIR}/${SAMPLE_ID}_1.fastq" \\
            -2 "${SAMPLE_DIR}/${SAMPLE_ID}_2.fastq" \\
            -S "${SAMPLE_DIR}/${SAMPLE_ID}_output.sam"
else
    # Single-end mapping
    bowtie2 -x "${SAMPLE_DIR}/${GENOME}_index" \\
            -U "${SAMPLE_DIR}/${SAMPLE_ID}.fastq" \\
            -S "${SAMPLE_DIR}/${SAMPLE_ID}_output.sam"
fi

# Convert SAM to BAM, sort it and remove intermediate files
samtools view -Sb "${SAMPLE_DIR}/${SAMPLE_ID}_output.sam" | samtools sort -o "${SAMPLE_DIR}/${SAMPLE_ID}_sorted_output.bam"
samtools index "${SAMPLE_DIR}/${SAMPLE_ID}_sorted_output.bam"

# Cleanup intermediate files
rm "${SAMPLE_DIR}/${SAMPLE_ID}_output.sam"

# Ensure coverm_output directory exists
mkdir -p "${SAMPLE_DIR}/coverm_output"

# Run inStrain for mapping profiles
inStrain quick_profile "$SAMPLE_DIR/${SAMPLE_ID}_sorted_output.bam" "$GENOME_PATH" \\
               -p 24 \\
               --breadth_cutoff 0.5 \\
               -o "$SAMPLE_DIR/${SAMPLE_ID}_instrain_output"

# Run coverm to calculate relative abundance
coverm genome -b "$SAMPLE_DIR/${SAMPLE_ID}_sorted_output.bam" \\
              -m relative_abundance \\
              --genome-fasta-directory "$SAMPLE_DIR" \\
              -t 24 \\
              --output-file "$SAMPLE_DIR/coverm_output/coverm_output.txt"

# Cleanup all other intermediate files
rm "${SAMPLE_DIR}/${SAMPLE_ID}_sorted_output.bam"
rm "${SAMPLE_DIR}/${SAMPLE_ID}_sorted_output.bam.bai"
# rm "${SAMPLE_DIR}/${GENOME}_index.*.bt2"

EOF
done
