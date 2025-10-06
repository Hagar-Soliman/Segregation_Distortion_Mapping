#!/bin/bash

# Create or overwrite the output job list file
output_file="step10_genotypes_list.txt"
> "$output_file"

# Loop through all matching files
for file in g.Genotypes.*.bam.txt; do
    # Extract sample ID by removing prefix and suffix
    sampleID=$(basename "$file" | sed 's/^g\.Genotypes\.//' | sed 's/\.bam\.txt$//')

    # Write the command line to the job list
    echo "python step10_calc_genotype_err_rate.py $sampleID" >> "$output_file"
done