#!/bin/bash
#SBATCH --job-name=UMB_dedup_P2_L1_py3
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=30G
#SBATCH --mail-user=piafranziska.schwarz@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#Script to remove PCR duplicates, step1 of processing ddRAD data for QTL mapping

module load miniconda
conda activate py3_env


INPUT_DIR=/home/ps2267/ycga_work/hybrid_necrosis/raw_data/Sample_PS_P2_L1_2024
OUTPUT_DIR=/home/ps2267/ycga_work/hybrid_necrosis/raw_data/Sample_PS_P2_L1_2024/rmdup
SCRIPT_DIR=/home/ps2267/ycga_work/hybrid_necrosis/raw_data/Sample_PS_P2_L1_2024

cd $INPUT_DIR

#unzip gzipped fastq files for python3 script
gzip -d PS_P*

#rename files from sequencing center to match format of R1, R2, i5 and i7 fastq files: 
#files from YCGA will be supplied as: 
#I1 = I7 index
#R1 = Sequencing Read 1
#R2 = I5 index read
#R3 = Sequencing Read 2

# rename library files to match format of R1, R2, i5 and i7 fastq files 

mv PS_P2_L1_2024_S1_L006_I1_001.fastq PS_P2_L1_2024_S1_L006_i7.fastq
mv PS_P2_L1_2024_S1_L006_R1_001.fastq PS_P2_L1_2024_S1_L006_R1.fastq
mv PS_P2_L1_2024_S1_L006_R2_001.fastq PS_P2_L1_2024_S1_L006_i5.fastq
mv PS_P2_L1_2024_S1_L006_R3_001.fastq PS_P2_L1_2024_S1_L006_R2.fastq


mkdir -p $OUTPUT_DIR

cd $INPUT_DIR

# run python3 script to remove duplicates, customize -p and -s arguments for each library 
python3 $SCRIPT_DIR/rmdup_molbarcodes_lila_r.py \
  -p PS_P2_L1_2024_S1_L006 \
  -s fastq

#gzip fastq files 
gzip $INPUT_DIR/*.rmdup_*.fastq


# check that files were generated
echo "Generated rmdup files:"
ls $INPUT_DIR/*.rmdup_*.fastq.gz

#move files to output directory
mv $INPUT_DIR/*.rmdup_*.fastq.gz $OUTPUT_DIR


conda deactivate
