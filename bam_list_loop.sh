#!/bin/bash

BASE_DIR="/home/hks25/palmer_scratch/libs/aligned"
PREFIX="1A"
OUTPUT_DIR="/home/hks25/palmer_scratch/SNP_calling"
OUTPUT_FILE="${OUTPUT_DIR}/${PREFIX}_bam_list.txt"

# Clear or create the output file
> "$OUTPUT_FILE"

# Loop through each ${PREFIX}*_aligned directory
for dir in "$BASE_DIR"/${PREFIX}*_aligned; do
    if [[ -d "$dir" ]]; then
        echo "Scanning: $dir"
        find "$dir" -type f -name "*.bam" >> "$OUTPUT_FILE"
    fi
done

echo "BAM list created at: $OUTPUT_FILE"