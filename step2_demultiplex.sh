#!/bin/bash
#SBATCH --job-name=1A_1A_demultiplex
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o demup_1A_1A.out
#SBATCH -e demup_1A_1A.err

###IMPORTANT:###
#before running this script, rename fastq.gz files to name format as received from the sequencing center. Otherwise stacks will not recognize the input files. 
#rename to format PS_P4_L1_2024_S3_L006.R1_001.fastq.gz (same for R2 file) or PS_P4_L1_2024_S3_L006_R1_001.fastq.gz (or with underscore before R#)

#note about the barcodes file: stacks is picky about the length of the file names for each sample. picky how many "columns" it can have (meeaning how many 
# sets of characters separated by _ (did not test if the issue is the overall number of characters or if _ creates a new "column" to the file. 

module load Stacks/2.59-GCCcore-10.2.0

library=1A_1A  #modify sample library name and make sure the paths are correct.
 
INPUT_DIR=/home/hks25/palmer_scratch/libs/${library}/rmdup
OUT_DIR=/home/hks25/palmer_scratch/libs/demultiplexed/${library}_demultiplexed
#barcodes=/home/hks25/ycga_work/hybrid_necrosis/demultiplexed/ #add file name to path at -b flag

mkdir -p /home/hks25/palmer_scratch/libs/demultiplexed/${library}_demultiplexed

# List files in the input directory
echo "Listing files in ${INPUT_DIR}:"
ls -l ${INPUT_DIR}

process_radtags -P -p ${INPUT_DIR} \
                -o ${OUT_DIR} \
                -b /home/hks25/palmer_scratch/libs/${library}/rmdup/${library}_barcodes.txt \
                --renz-1 PstI --renz-2 BfaI -r -c -q --inline_null --bestrad \
                --rescue


