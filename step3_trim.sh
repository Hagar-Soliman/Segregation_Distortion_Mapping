#!/bin/bash
#SBATCH --job-name=1A_1A_trimmed
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=16G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


#run this script following demultiplexing with stacks

module load Trimmomatic/0.39-Java-11

library=1A_1A
echo "library: $library"

input_dir=/home/hks25/palmer_scratch/libs/demultiplexed/${library}_demultiplexed #directory that contains files from demultiplexing with stacks 
output_dir=/home/hks25/palmer_scratch/libs/trimmed/${library}_trimmed #directory to store trimmed files

mkdir -p $output_dir

#Search all .1.fq.gz files without .rem in the demultiplexed_output directory and look for the corresponding .2.fq.gz file.

for file_R1 in $(ls $input_dir/*_*.1.fq.gz | grep -v "\.rem\.") #grabs files in the format of stacks output 
do
    base_name=$(basename $file_R1 ".1.fq.gz")
    
    file_R2="$input_dir/${base_name}.2.fq.gz"

    input_file_R1="$input_dir/${base_name}.1.fq.gz"
    input_file_R2="$input_dir/${base_name}.2.fq.gz"

    output_file_R1_paired="$output_dir/${base_name}.1_paired_trimmed.fq.gz"
    output_file_R1_unpaired="$output_dir/${base_name}.1_unpaired_trimmed.fq.gz"
    output_file_R2_paired="$output_dir/${base_name}.2_paired_trimmed.fq.gz"
    output_file_R2_unpaired="$output_dir/${base_name}.2_unpaired_trimmed.fq.gz"

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 16 \
        $input_file_R1 $input_file_R2 \
        $output_file_R1_paired $output_file_R1_unpaired \
        $output_file_R2_paired $output_file_R2_unpaired \
        ILLUMINACLIP:/home/hks25/palmer_scratch/step3_adapters.fa:2:30:10:2:TRUE \
        LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:30
done
