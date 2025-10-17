#!/bin/bash

## more quality control post quality checks
# SLURM settings for job creation
FASTQ_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/bbduk_quality_processed/human_reads_removed"
RESULTS_DIR="./read_files_putative_mucins/testing_if_we_should_subsample/fastqc_results"
job_dir="./read_files_putative_mucins/testing_if_we_should_subsample/fastqc_scripts"

# Create the results and job directories if they don't exist
mkdir -p "$RESULTS_DIR"
mkdir -p "$job_dir"

# Function to create a unique SLURM script for each file pair
create_fastqc_script() {
    fastq_dir=$1
    fastq_id=$(basename "$fastq_dir")
    fastq_file_1="${fastq_dir}/${fastq_id}_1.fastq"
    fastq_file_2="${fastq_dir}/${fastq_id}_2.fastq"
    job_script="${job_dir}/fastqc_${fastq_id}.sh"

    cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH --time=30:00:00           # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=70G        # Maximum memory required per CPU
#SBATCH --job-name=fastqc_${fastq_id}
#SBATCH --output=${RESULTS_DIR}/fastqc_${fastq_id}.%J.out
#SBATCH --error=${RESULTS_DIR}/fastqc_${fastq_id}.%J.err
#SBATCH --partition=yinlab,batch,guest

ml fastqc/0.12  # Load FastQC

fastq_file_1="${fastq_file_1}"
fastq_file_2="${fastq_file_2}"

echo "Running FastQC on \$fastq_file_1 and \$fastq_file_2..."

# Run FastQC on both paired files
fastqc "\$fastq_file_1" "\$fastq_file_2" --outdir="${RESULTS_DIR}"

# Unzip the FastQC results
unzip "${RESULTS_DIR}/${fastq_id}_1_fastqc.zip" -d "${RESULTS_DIR}"
unzip "${RESULTS_DIR}/${fastq_id}_2_fastqc.zip" -d "${RESULTS_DIR}"

# Check if "Per base sequence quality" and "Per sequence quality scores" passed for both files
if grep -q "FAIL" "${RESULTS_DIR}/${fastq_id}_1_fastqc/summary.txt" | grep -E "Per base sequence quality|Per sequence quality scores" || grep -q "FAIL" "${RESULTS_DIR}/${fastq_id}_2_fastqc/summary.txt" | grep -E "Per base sequence quality|Per sequence quality scores"; then
    echo "One or both of the paired files (${fastq_id}_1.fastq and ${fastq_id}_2.fastq) failed either 'Per base sequence quality' or 'Per sequence quality scores'."
else
    echo "Both paired files (${fastq_id}_1.fastq and ${fastq_id}_2.fastq) passed 'Per base sequence quality' and 'Per sequence quality scores'."
fi

# Clean up FastQC files after logging the results
rm -r "${RESULTS_DIR}/${fastq_id}_1_fastqc" "${RESULTS_DIR}/${fastq_id}_2_fastqc"
rm "${RESULTS_DIR}/${fastq_id}_1_fastqc.zip" "${RESULTS_DIR}/${fastq_id}_2_fastqc.zip"
EOF
}

# Find all directories containing FASTQ files and generate SLURM scripts for them
find "$FASTQ_DIR" -mindepth 1 -maxdepth 1 -type d | while read fastq_dir; do
    create_fastqc_script "$fastq_dir"
done

echo "FastQC SLURM scripts created in $job_dir."
