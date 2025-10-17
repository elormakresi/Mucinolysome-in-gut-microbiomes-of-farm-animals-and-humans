#!/bin/bash

## to subsample the reads for space/memory

# SLURM settings for job creation
input_dir="./read_files_putative_mucins/reads_for_prevalence/extract_from_attic/jerry_swan_work/canine_gut/bbduk_quality_processed/human_reads_removed"
output_base_dir="./read_files_putative_mucins/reads_for_prevalence/extract_from_attic/jerry_swan_work/canine_gut/bbduk_quality_processed/human_reads_removed/subsampled_5mil"
target_reads=10000000
seed=42
job_dir="./read_files_putative_mucins/reads_for_prevalence/extract_from_attic/jerry_swan_work/subsample_slurm_jobs_caninegut"

# Create a directory for the SLURM job scripts
mkdir -p "$job_dir"
mkdir -p "$output_base_dir"

# Function to create a unique SLURM script for each file pair
create_slurm_script() {
    fastq_file_1=$1
    base_name=$(basename "$fastq_file_1" _1.fastq)
    fastq_file_2="${base_name}_2.fastq"

    job_script="${job_dir}/subsample_${base_name}.sh"

    # Generate a unique SLURM script for the current FASTQ file pair
    cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH --time=40:00:00           # Run time in hh:mm:ss
#SBATCH --mem=600G                # Memory required
#SBATCH --cpus-per-task=2         # Request 1 core
#SBATCH --job-name=subsample_${base_name}
#SBATCH --output=subsample_${base_name}.%J.out
#SBATCH --error=subsample_${base_name}.%J.err
#SBATCH --partition=yinlab,batch,guest

# Source the bashrc file to ensure conda is initialized
source ~/.bashrc

# Activate the conda environment
conda activate subsample

# Set the path to the paired FASTQ files
fastq_file_1="${input_dir}/${base_name}_1.fastq"
fastq_file_2="${input_dir}/${base_name}_2.fastq"

# Automatically create the output directory
output_dir="${output_base_dir}/${base_name}"
mkdir -p "\$output_dir"

# Output files
output_file_1="\$output_dir/${base_name}_1.fastq"
output_file_2="\$output_dir/${base_name}_2.fastq"

# Check the number of reads and subsample if needed
read_count=\$(seqkit stats -T "\$fastq_file_1" | awk 'NR==2 {print \$4}')
if (( read_count >= $target_reads )); then
    echo "Subsampling \$fastq_file_1 and \$fastq_file_2 to $target_reads reads..."

    # Subsample with seqkit
    seqkit sample -n $target_reads -s $seed "\$fastq_file_1" -o "\$output_file_1"
    seqkit sample -n $target_reads -s $seed "\$fastq_file_2" -o "\$output_file_2"
else
    echo "Skipping \$fastq_file_1 and \$fastq_file_2 because they have fewer than $target_reads reads."
fi

echo "Subsampling complete for ${base_name}"
EOF
}

# Find all *_1.fastq files and generate SLURM scripts for them
find "$input_dir" -name "*_1.fastq" | while read fastq_file_1; do
    create_slurm_script "$fastq_file_1"
done

echo "SLURM scripts created in $job_dir"
