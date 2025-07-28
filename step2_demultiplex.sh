#!/bin/bash
#SBATCH --job-name=p7l2_demultiplex
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --mail-user=piafranziska.schwarz@yale.edu
#SBATCH -o demu_p7l2.out
#SBATCH -e demu_p7l2.err

###IMPORTANT:###
#before running this script, rename fastq.gz files to name format as received from the sequencing center. Otherwise stacks will not recognize the input files. 
#rename to format PS_P4_L1_2024_S3_L006.R1_001.fastq.gz (same for R2 file) or PS_P4_L1_2024_S3_L006_R1_001.fastq.gz (or with underscore before R#)

#note about the barcodes file: stacks is picky about the length of the file names for each sample. picky how many "columns" it can have (meeaning how many 
# sets of characters separated by _ (did not test if the issue is the overall number of characters or if _ creates a new "column" to the file. 

module load Stacks/2.59-GCCcore-10.2.0

library=Sample_PS_P7_L2_2024  #modify sample library name
 
MASTER_DIR=/home/ps2267/ycga_work/hybrid_necrosis/raw_data/${library}/
INPUT_DIR=/home/ps2267/ycga_work/hybrid_necrosis/raw_data/${library}/rmdup
OUT_DIR=/home/ps2267/ycga_work/hybrid_necrosis/demultiplexed/${library}_demultiplexed
#barcodes=/home/ps2267/ycga_work/hybrid_necrosis/demultiplexed/ #add file name to path at -b flag

#mkdir -p /home/ps2267/ycga_work/hybrid_necrosis/demultiplexed
mkdir /home/ps2267/ycga_work/hybrid_necrosis/demultiplexed/${library}_demultiplexed

# List files in the input directory
echo "Listing files in ${INPUT_DIR}:"
ls -l ${INPUT_DIR}

#cd ${MASTER_DIR}/

process_radtags -P -p ${INPUT_DIR} \
                -o ${OUT_DIR} \
                -b /home/ps2267/ycga_work/hybrid_necrosis/demultiplexed/P7_L2_sample_barcodes.txt \
                --renz-1 PstI --renz-2 BfaI -r -c -q --inline_null --bestrad \
                --rescue


