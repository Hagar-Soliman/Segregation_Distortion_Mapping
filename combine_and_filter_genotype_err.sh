#!/bin/bash

# Combines step10 output files and filters individuals with e1, e2, beta all < 0.2
# Output: tab-delimited file with header

COMBINED="5C_combined_genotypes_err.txt"
FILTERED="5C_filtered_genotypes_err.txt"

# --- Step 1: Combine ---
> "$COMBINED"
for file in dsq*.out; do
    if [[ -f "$file" && -s "$file" ]]; then
        grep '^5C_' "$file" | tr -s ' ' '\t' >> "$COMBINED"
    fi
done

echo "Combined: $(wc -l < "$COMBINED") individuals"

# --- Step 2: Filter (e1 < 0.2, e2 < 0.2, beta < 0.2) ---
awk -F'\t' '$3 < 0.2 && $4 < 0.2 && $5 < 0.2' "$COMBINED" > "$FILTERED"

echo "Filtered: $(( $(wc -l < "$FILTERED") - 1 )) individuals passed (e1 < 0.2, e2 < 0.2, beta < 0.2)"
echo "Dropped:  $(( $(wc -l < "$COMBINED") - $(wc -l < "$FILTERED") + 1 )) individuals"

# --- Step 3: Extract retained individual names ---
INDIV_LIST="5C_genotype_err_indivi_list.txt"
awk -F'\t' '{print $1}' "$FILTERED" > "$INDIV_LIST"
echo "Individual list written: $INDIV_LIST ($(wc -l < "$INDIV_LIST") names)"
