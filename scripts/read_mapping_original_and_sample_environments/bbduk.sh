#!/bin/bash

## to perform  quality control. Inputs are fastq files 

# SLURM settings for job creation
PARENT_DIR="./read_files_putative_mucins/testing_if_we_should_subsample"
BBDUK_PATH="./bbmap/bbduk.sh"
ADAPTERS_FA="./bbmap/resources/adapters.fa"
job_dir="./read_files_putative_mucins/testing_if_we_should_subsample/bbduk_slurm"

# Create the job directory if it doesn't exist
mkdir -p "$job_dir"

# Function to create a unique SLURM script for each sample
create_bbduk_script() {
    sample_dir=$1
    sample_name=$(basename "$sample_dir")
    job_script="${job_dir}/bbduk_${sample_name}.sh"

    input_R1="${sample_dir%/}/${sample_name}_1.fastq"
    input_R2="${sample_dir%/}/${sample_name}_2.fastq"
    output_R1="${sample_dir%/}/${sample_name}_trimmed_1.fastq"
    output_R2="${sample_dir%/}/${sample_name}_trimmed_2.fastq"

    # Generate the SLURM script for BBDuk processing
    cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH --time=80:00:00           # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=200G         # Maximum memory required per CPU
#SBATCH --job-name=hunterbbduk_${sample_name}
#SBATCH --output=hunterbbduk_${sample_name}.%J.out
#SBATCH --error=hunterbbduk_${sample_name}.%J.err
#SBATCH --partition=yinlab,batch,guest

ml java/12  # Load Java

# Run BBDuk for adapter trimming, quality trimming, and length filtering
"$BBDUK_PATH" in1="$input_R1" in2="$input_R2" out1="$output_R1" out2="$output_R2" \\
ref="$ADAPTERS_FA" ktrim=r k=23 mink=11 hdist=1 \\
trimq=16 minlen=55 tbo tpe

echo "Processing completed for $sample_name"
EOF
}

# Loop through each sample directory and create SLURM scripts
for sample_dir in "$PARENT_DIR"/*/; do
    create_bbduk_script "$sample_dir"
done

echo "BBDuk SLURM scripts created in $job_dir."
