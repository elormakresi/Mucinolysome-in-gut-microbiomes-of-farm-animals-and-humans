#!/bin/bash
#SBATCH --time=15:00:00      
#SBATCH --mem=20gb
#SBATCH --job-name=trimmomatic
#SBATCH --error=trimmomatic.err
#SBATCH --output=trimmomatic.out

module load java 
module load trimmomatic/0.39 

cd ./01.RawData 

INPATH=./01.RawData/trim

trimmomatic PE G1/G1_1.fq G1/G1_2.fq $INPATH/G1_1_p.fq $INPATH/G1_1_up.fq $INPATH/G1_2_p.fq $INPATH/G1_2_up.fq \
ILLUMINACLIP:/util/opt/anaconda/deployed-conda-envs/packages/trimmomatic/envs/trimmomatic-0.39/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:126 CROP:125

trimmomatic PE G2/G2_1.fq G2/G2_2.fq $INPATH/G2_1_p.fq $INPATH/G2_1_up.fq $INPATH/G2_2_p.fq $INPATH/G2_2_up.fq \
ILLUMINACLIP:/util/opt/anaconda/deployed-conda-envs/packages/trimmomatic/envs/trimmomatic-0.39/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:126 CROP:125

trimmomatic PE G3/G3_1.fq G3/G3_2.fq $INPATH/G3_1_p.fq $INPATH/G3_1_up.fq $INPATH/G3_2_p.fq $INPATH/G3_2_up.fq \
ILLUMINACLIP:/util/opt/anaconda/deployed-conda-envs/packages/trimmomatic/envs/trimmomatic-0.39/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:126 CROP:125

trimmomatic PE M1/M1_1.fq M1/M1_2.fq $INPATH/M1_1_p.fq $INPATH/M1_1_up.fq $INPATH/M1_2_p.fq $INPATH/M1_2_up.fq \
ILLUMINACLIP:/util/opt/anaconda/deployed-conda-envs/packages/trimmomatic/envs/trimmomatic-0.39/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:126 CROP:125

trimmomatic PE M2/M2_1.fq M2/M2_2.fq $INPATH/M2_1_p.fq $INPATH/M2_1_up.fq $INPATH/M2_2_p.fq $INPATH/M2_2_up.fq \
ILLUMINACLIP:/util/opt/anaconda/deployed-conda-envs/packages/trimmomatic/envs/trimmomatic-0.39/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:126 CROP:125

trimmomatic PE M3/M3_1.fq M3/M3_2.fq $INPATH/M3_1_p.fq $INPATH/M3_1_up.fq $INPATH/M3_2_p.fq $INPATH/M3_2_up.fq \
ILLUMINACLIP:/util/opt/anaconda/deployed-conda-envs/packages/trimmomatic/envs/trimmomatic-0.39/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:126 CROP:125
