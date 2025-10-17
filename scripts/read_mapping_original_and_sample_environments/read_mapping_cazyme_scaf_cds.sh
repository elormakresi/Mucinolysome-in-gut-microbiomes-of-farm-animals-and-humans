#!/bin/bash

# Directories
MAG_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/scaffoldin_CDS_all"
SAMPLES_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/bbduk_quality_processed/human_reads_removed/high_quality"
SLURM_SCRIPTS_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/read_mapping_scripts_scaf"

# Create the SLURM scripts directory
mkdir -p "$SLURM_SCRIPTS_DIR"

for MAG_FOLDER in "$MAG_DIR"/*; do
    MAG_ID=$(basename "$MAG_FOLDER")
    MAG="$MAG_FOLDER/${MAG_ID}.fna"

    for SAMPLE_DIR in "$SAMPLES_DIR"/*; do
        SAMPLE_ID=$(basename "$SAMPLE_DIR")
        OUTPUT_DIR="$MAG_FOLDER/$SAMPLE_ID"
        mkdir -p "$OUTPUT_DIR/coverm_output"  # Create the coverm_output directory

        # Define SLURM script content
        cat <<EOF > "$SLURM_SCRIPTS_DIR/${MAG_ID}_${SAMPLE_ID}.sh"
#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mem=120gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=${MAG_ID}_${SAMPLE_ID}
#SBATCH --partition=yinlab,batch,guest

ml bowtie/2.3
ml instrain/1.3
ml coverm/0.6
ml samtools/1.19

# Build Bowtie2 index
bowtie2-build "$MAG" "$OUTPUT_DIR/${MAG_ID}_index"

# added --very-sensitive cos default had low alignment rates removed on may_9th_2025
if [ -f "${SAMPLE_DIR}/${SAMPLE_ID}_2.fastq" ]; then
    # Paired-end mapping
    bowtie2 -x "$OUTPUT_DIR/${MAG_ID}_index" \
            -1 "${SAMPLE_DIR}/${SAMPLE_ID}_1.fastq" \
            -2 "${SAMPLE_DIR}/${SAMPLE_ID}_2.fastq" \
            -p 5 \
            -S "$OUTPUT_DIR/${SAMPLE_ID}_output.sam"
else
    # Single-end mapping
    bowtie2 -x "$OUTPUT_DIR/${MAG_ID}_index" \
            -U "${SAMPLE_DIR}/${SAMPLE_ID}.fastq" \
            -p 5 \
            -S "$OUTPUT_DIR/${SAMPLE_ID}_output.sam"
fi

# Convert SAM to BAM, sort, and index
samtools view -Sb "$OUTPUT_DIR/${SAMPLE_ID}_output.sam" > "$OUTPUT_DIR/${SAMPLE_ID}_output.bam"
samtools sort "$OUTPUT_DIR/${SAMPLE_ID}_output.bam" -o "$OUTPUT_DIR/${SAMPLE_ID}_sorted_output.bam"
samtools index "$OUTPUT_DIR/${SAMPLE_ID}_sorted_output.bam"

# Remove SAM to save space
rm "$OUTPUT_DIR/${SAMPLE_ID}_output.sam"

# Run inStrain quick_profile for breadth & coverage
inStrain quick_profile "$OUTPUT_DIR/${SAMPLE_ID}_sorted_output.bam" "$MAG" \
               -p 5 \
               --breadth_cutoff 0.5 \
               -o "$OUTPUT_DIR/${SAMPLE_ID}_instrain_output"

# Run inStrain profile for detailed analysis
inStrain profile "$OUTPUT_DIR/${SAMPLE_ID}_sorted_output.bam" "$MAG" \
               -p 5 \
               -o "$OUTPUT_DIR/${SAMPLE_ID}_instrain_output"

# Run coverm for relative abundance
coverm genome -b "$OUTPUT_DIR/${SAMPLE_ID}_sorted_output.bam" \
              --min-covered-fraction 0.5 \
              -m relative_abundance \
              --genome-fasta-directory "$MAG_FOLDER" \
              -t 5 \
              --output-file "$OUTPUT_DIR/coverm_output/coverm_output.txt"

EOF
    done
done
