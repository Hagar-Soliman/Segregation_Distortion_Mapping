#!/bin/bash

library="1A"
# Define base path
base_dir="/home/hks25/palmer_scratch/libs/aligned"

# Define destination directory
dest_dir="$base_dir/${library}_refalt"

# Create destination directory if it doesn't exist
mkdir -p "$dest_dir"

# Loop through all subdirectories starting with 1A_
for dir in "$base_dir"/${library}_*/; do
  # Find and move all .bam.txt files in each subdirectory
  find "$dir" -maxdepth 1 -type f -name "*.bam.txt" -exec mv {} "$dest_dir" \;
done

echo "All .bam.txt files moved to $dest_dir"