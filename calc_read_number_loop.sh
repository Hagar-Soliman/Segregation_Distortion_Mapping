#!/bin/bash
#SBATCH --job-name=calc_read_number_1A1A
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err





module load SAMtools

echo "" > all.stats

for file in *bam ; do
reads=$(samtools flagstat $file | head -1)
echo "$file $reads" >> all.stats

done