#!/bin/bash
#SBATCH --job-name=rmdup_1A_1A
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=30G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#Script to remove PCR duplicates, step1 of processing ddRAD data for mapping

module load miniconda
conda activate py3_env


INPUT_DIR=/home/hks25/palmer_scratch/libs/1A_1A/Unaligned
OUTPUT_DIR=/home/hks25/palmer_scratch/libs/1A_1A/rmdup
SCRIPT_DIR=/home/hks25/palmer_scratch

cd $INPUT_DIR

#unzip gzipped fastq files for python3 script
gzip -d 1A_1A*

#rename files from sequencing center to match format of R1, R2, i5 and i7 fastq files: 
#files from YCGA will be supplied as: 
#I1 = i7 index
#R1 = Sequencing Read 1
#R2 = i5 index read
#R3 = Sequencing Read 2

# rename library files to match format of R1, R2, i5 and i7 fastq files 

mv 1A_1A_S1_L003_I1_001.fastq 1A_1A_i7.fastq
mv 1A_1A_A22VH7TLT4_L003_R1_001.fastq 1A_1A_R1.fastq
mv 1A_1A_A22VH7TLT4_L003_R2_001.fastq 1A_1A_i5.fastq
mv 1A_1A_A22VH7TLT4_L003_R3_001.fastq 1A_1A_R2.fastq


mkdir -p $OUTPUT_DIR

cd $INPUT_DIR

# run python3 script to remove duplicates, customize -p and -s arguments for each library 
python3 $SCRIPT_DIR/step1_rmdup.py \
  -p 1A_1A \
  -s fastq

#gzip fastq files 
gzip $INPUT_DIR/*.rmdup_*.fastq


# check that files were generated
echo "Generated rmdup files:"
ls $INPUT_DIR/*.rmdup_*.fastq.gz

#move files to output directory
mv $INPUT_DIR/*.rmdup_*.fastq.gz $OUTPUT_DIR


conda deactivate
