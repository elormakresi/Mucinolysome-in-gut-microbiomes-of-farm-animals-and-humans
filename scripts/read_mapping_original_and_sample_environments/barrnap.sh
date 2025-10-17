#!/bin/bash
#SBATCH --time=5:00:00           # Run time in hh:mm:ss
#SBATCH --mem=60G                # Memory required
#SBATCH --cpus-per-task=2         # Request 2 cores
#SBATCH --job-name=barrnap_default
#SBATCH --output=barrnap.%J.out
#SBATCH --error=barrnap.%J.err
#SBATCH --partition=yinlab,batch,guest


## to predict 16s RNAs in the 65 genomes
# Define input and output paths
INPATH=./all_fna_seqs_uniq
OUTPATH=./16s_seqs_final
SEQ_OUTPATH=$OUTPATH/extracted_16s

# Ensure output directories exist
mkdir -p $OUTPATH
mkdir -p $SEQ_OUTPATH

# Activate Barrnap environment
source ~/.bashrc
conda activate barrnap_env

# Loop through each .fna file in the input directory
for mag in $INPATH/*.fna; do
    # Extract the base name of the file (without the full path)
    base=$(basename $mag .fna)

    # Run Barrnap with stricter defaults for MAGs
    barrnap --kingdom bac \
        --outseq $SEQ_OUTPATH/${base}_16s.fna \
        $mag > $OUTPATH/${base}_all.gff

    # Filter to retain only 16S rRNA entries in the GFF
    grep "16S_rRNA" $OUTPATH/${base}_all.gff > $OUTPATH/${base}_16s.gff
done

# Deactivate Conda environment after the loop
conda deactivate
