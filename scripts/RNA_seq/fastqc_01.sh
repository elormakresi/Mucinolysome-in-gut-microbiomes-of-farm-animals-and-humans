#!/bin/bash
#SBATCH --time=01:00:00      
#SBATCH --mem=10gb
#SBATCH --job-name=fastqc
#SBATCH --error=fastqc.err
#SBATCH --output=fastqc.out

ml fastqc/0.12
ml multiqc/py37/1.8

fastqc -o ./01.RawData/fastqc_output ./01.RawData/G1/G1_1.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/G1/G1_2.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/G2/G2_1.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/G2/G2_2.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/G3/G3_1.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/G3/G3_2.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/M1/M1_1.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/M1/M1_2.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/M2/M2_1.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/M2/M2_2.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/M3/M3_1.fq
fastqc -o ./01.RawData/fastqc_output ./01.RawData/M3/M3_2.fq

multiqc ./01.RawData/fastqc_output 

