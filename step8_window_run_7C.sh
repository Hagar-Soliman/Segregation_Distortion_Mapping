#!/bin/bash
#SBATCH --job-name=7C_step8_windows
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

map="7C"
module load Python/3.12.3-GCCcore-13.3.0

input_directory="/home/hks25/palmer_scratch/googaV3/refAlt/${map}_refAlt"
script_directory="/home/hks25/palmer_scratch/googaV3/refAlt/${map}_refAlt"

for window_size in 50kb; do
    output_directory="${input_directory}/${window_size}_windows"
    mkdir -p "$output_directory"

    for file in "$input_directory"/*.bam.txt; do
        base_name=$(basename "$file")

        if [[ "$base_name" == "IM62.bam.txt" || "$base_name" == "CCC9_GDS.bam.txt" ]]; then
            continue
        fi

        name_without_txt="${base_name%.txt}"
        echo "Processing $file with ${window_size} windows"
        python3 "${script_directory}/step8_50kb_window.py" "$file" "$name_without_txt" "$output_directory"
    done
done
