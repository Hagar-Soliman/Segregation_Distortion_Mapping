#!/bin/bash
#SBATCH --job-name=1A_step8_windows
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

library="1A"
module load Python/3.12.3-GCCcore-13.3.0

input_directory="/home/hks25/palmer_scratch/libs/aligned/${library}_refalt"
output_directory="/home/hks25/palmer_scratch/libs/aligned/${library}_refalt/windowpy"

mkdir -p "$output_directory"

for file in "$input_directory"/*.bam.txt
do

  if [ -f "$file" ]; then
    base_name=$(basename "$file")
    name_without_txt="${base_name%.txt}"
    echo "Processing file: $file with analysis name: $name_without_txt"

    python3 /home/hks25/palmer_scratch/step8_window.py "$file" "$name_without_txt" "$output_directory"
  fi
done
