#!/bin/bash
#SBATCH --job-name=1A_1A_align
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL 
#SBATCH --mem=64G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



module load BWA/0.7.17-GCCcore-12.2.0
module load SAMtools/1.20-GCC-12.2.0

#Change library name and make sure paths are correct!!!!!!!!!!!!!!!!!!11

library=1A_1A

echo "library: $library"

input_dir=/home/hks25/palmer_scratch/libs/trimmed/${library}_trimmed
output_dir=/home/hks25/palmer_scratch/libs/aligned/${library}_aligned

mkdir -p "$output_dir"

REF_PATH=/gpfs/gibbs/project/coughlan/shared/genomes/ref/Mimulus_guttatus_var_IM62_v3.mainGenome.fasta

for file_1 in $input_dir/*.1_paired_trimmed.fq.gz
do
    # Get the base name of the file (excluding path and suffix)
    base_name=$(basename $file_1 ".1_paired_trimmed.fq.gz")
    file_2="${input_dir}/${base_name}.2_paired_trimmed.fq.gz"

    # checks if the paired file file_2 exists. The -f option tests if the file is a regular file.
    if [[ -f $file_2 ]]; then
        # following line defines the output SAM file path
        output_bam="${output_dir}/${base_name}.bam"

        # following line runs BWA alignement using 16 threads and outputs SAM files to output directory
       bwa mem -t 16 $REF_PATH $file_1 $file_2 | samtools view -@ 16 -b -F 4 -q 20 | samtools sort -@ 16 -o "$output_bam"
       
        # Index the BAM file after it's created
        samtools index "$output_bam"
    else
        echo "Pairing-file ${file_2} does not exist, skip ${file_1}"
    fi
done
