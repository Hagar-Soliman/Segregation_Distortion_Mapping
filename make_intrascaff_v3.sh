#!/bin/bash

# =============================================================================
# make_intrascaff_v3.sh
#
# PURPOSE:
#   Converts the intrascaff v2 file (recombination fractions per interval)
#   into the v3 format required by step 13, by adding a 4th column containing
#   cumulative genetic map positions in cM (centiMorgans).
#
#   The conversion follows the same logic as the manual Excel method:
#     - First marker of each chromosome:  cM = 0
#     - Each subsequent marker j:         cM = cM[j-1] + 100 * r[j-1]
#       where r[j-1] is the recombination fraction of the previous marker
#     - The last marker per chromosome has a large negative value (log-likelihood)
#       in column 3 — this is not a recombination rate and is NOT added to cM.
#       Its cM value just carries forward from the previous marker.
#
#   The output preserves all original columns and appends the cM as column 4.
#   Output is tab-delimited with no header, as required by step 13.
#
# USAGE:
#   bash make_intrascaff_v3.sh
#   (auto-detects v2 file matching *intrascaff*v2*.txt in current directory)
#
# OUTPUT:
#   1A_intrascaff_v3.txt — tab-delimited, 4 columns:
#     chromosome  position  recomb_fraction  cumulative_cM
# =============================================================================

# Auto-detect v2 input file
V2=$(ls *intrascaff*v2*.txt 2>/dev/null | head -1)
if [[ -z "$V2" ]]; then
    echo "ERROR: No v2 intrascaff file found matching *intrascaff*v2*.txt"
    exit 1
fi

# Derive output filename by replacing v2 with v3
V3="${V2/v2/v3}"
# Strip any window-size prefix from the name to match expected step13 filename
# e.g. 1A_intrascaff50kb_v2.txt → 1A_intrascaff_v3.txt
V3=$(echo "$V3" | sed 's/intrascaff[0-9]*kb/intrascaff/')

echo "Input:  $V2"
echo "Output: $V3"
echo ""

# -----------------------------------------------------------------------------
# Core conversion using awk:
# - Track the current chromosome to detect when a new one starts (reset cM=0)
# - For each line, output the original 3 columns plus the running cM total
# - Only add r to the cumulative total if r >= 0 (skip the LL on the last line)
# -----------------------------------------------------------------------------

awk -F'\t' 'BEGIN { OFS="\t"; prev_scaff=""; cum_cM=0; prev_r=0 }
{
    scaff = $1
    pos   = $2
    r     = $3

    # New chromosome detected — reset cumulative cM and previous rate to 0
    if (scaff != prev_scaff) {
        cum_cM  = 0
        prev_r  = 0
        prev_scaff = scaff
    } else {
        # Add the previous marker recombination fraction to cumulative cM
        # Only if previous r >= 0 (negative value = log-likelihood on last marker)
        if (prev_r + 0 >= 0) {
            cum_cM += 100 * prev_r
        }
    }

    print scaff, pos, r, cum_cM

    prev_r = r

}' "$V2" > "${V3}.tmp"

# -----------------------------------------------------------------------------
# Sort chromosomes numerically by the number after "Chr_" (e.g. Chr_01...Chr_14)
# The cM computation above must happen first (it is order-dependent within each
# chromosome), so we sort only after the cumulative cM has been calculated.
# Within each chromosome, the original marker order is preserved by sort's
# stability (-s) combined with the secondary sort on position (-k2,2n).
# -----------------------------------------------------------------------------
sort -t$'\t' -k1,1 -k2,2n -s \
    --key=1,1 \
    -V \
    "${V3}.tmp" > "$V3"

rm "${V3}.tmp"

echo "Done. Lines written: $(wc -l < "$V3")"
echo ""
echo "First 5 lines:"
head -5 "$V3"
echo "..."
echo "Last 5 lines:"
tail -5 "$V3"
