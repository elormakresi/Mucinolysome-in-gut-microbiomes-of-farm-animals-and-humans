#!/usr/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=264gb
#SBATCH --ntasks=8
#SBATCH --job-name=eggnog
#SBATCH --partition=yinlab,batch,guest
#SBATCH --error=eggnog.err
#SBATCH --output=eggnog.out

ml eggnog-mapper/2.1

## functional annotation of kol109 reference 

# Copy the eggnog data directory to shared memory for faster access
cp -r $EGGNOG_DATA_DIR /dev/shm/eggnog_data_dir/

cd ./01.RawData/KO_annotations_ref_based

# Run eggnog-mapper
emapper.py \
  -i ./faa_files_putative_mucins/kol109.faa \
  -o kol109_eggnog \
  --itype proteins \
  --cpu 8 \
  --data_dir /dev/shm/eggnog_data_dir \
  --dbmem
