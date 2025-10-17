#!/usr/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=128gb
#SBATCH --ntasks=8
#SBATCH --job-name=read_mapping
#SBATCH --partition=yinlab,batch,guest
#SBATCH --error=feature_counts.err
#SBATCH --output=feature_counts.out


ml subread/2.0 

GTF_DIR="./reference_based_rna_seq"
BAM_DIR="./reference_based_rna_seq/output_bam_aligned"
OUTPUT_COUNTS="gene_counts.txt"

GTF_FILE="${GTF_DIR}/kol109.gtf"

echo "Collecting BAM files for featureCounts..."

BAM_FILES=$(ls "${BAM_DIR}"/*_sorted.bam)

# Convert the list to a space-separated string for featureCounts
BAM_FILES_STRING=""
for bam_file in $BAM_FILES; do
    BAM_FILES_STRING+=" ${bam_file}"
done

echo "Running featureCounts with the following BAM files:"
echo "$BAM_FILES_STRING"
echo ""

featureCounts -a "$GTF_FILE" \
              -F GTF \
              -t transcript \
              -g gene_id \
              -o "$OUTPUT_COUNTS" \
              -p \
              -s 2 \
              $BAM_FILES_STRING \
              --verbose

# Check if featureCounts was successful
if [ $? -ne 0 ]; then
    echo "featureCounts failed. Exiting."
    exit 1
fi

echo "Gene counts saved to ${OUTPUT_COUNTS}"
echo "FeatureCounts summary output (e.g., total reads, assigned reads) can be found in ${OUTPUT_COUNTS}.summary"
