#!/bin/bash
#SBATCH --job-name=1A_1A_sort_clean
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=32G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load SAMtools/1.20-GCC-12.2.0

library=1A_1A #change the library name and make sure the paths are correct
echo "library: $library"

input_dir=/home/hks25/palmer_scratch/libs/aligned_BWA/${library}_aligned
output_dir=/home/hks25/palmer_scratch/libs/sort_clean/${library}_sort_clean 

mkdir -p $output_dir

#search all SAM files, sort and clean them
for sam_file in $input_dir/*.sam
do
    #Determine the base name of the file (without path and suffix).
    base_name=$(basename $sam_file ".sam")
    
    # Output the sorted BAM file path
    output_bam="${output_dir}/${base_name}_sorted.bam"
    output_clean_bam="${output_dir}/${base_name}_sorted_clean.bam"

    #Convert SAM to BAM with Samtools and sort it
    samtools sort -@ 16 -o $output_bam $sam_file

    #Remove the original SAM file to conserve disk space
    rm $sam_file

    #Use Samtools for filtering: Remove unmapped reads and reads with an alignment quality below 30
    #apply -q 30 for downstream analysis of variant calling, could use -q 20 for exploratory analysis 
    #Pei-wei mentioned that 30 is too low and it removed a lot of reads.
    samtools view -@ 16 -b -F 4 -q 20 $output_bam > $output_clean_bam

    #Delete the unfiltered BAM to conserve disk space
    rm $output_bam
    
    echo "${base_name} sorting and cleaning have been completed, generated ${output_clean_bam}"
done
