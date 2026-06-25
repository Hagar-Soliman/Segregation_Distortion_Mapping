#!/bin/bash

# =============================================================================
# check_bad_markers.sh
#
# PURPOSE:
#   Two-round bad marker detection for the GOOGA pipeline.
#
#   ROUND 1 (run after step 11, on v1 intrascaff):
#     Catches markers where the HMM hit the 0.25 MaxRR ceiling — a hard sign
#     that one flanking marker is bad (paralogous, repeat, mapping artifact).
#
#   ROUND 2 (run after step 12, on v2 intrascaff):
#     Uses a lower threshold (RR_THRESHOLD = 0.1) to catch markers that were
#     not at the ceiling but still have suspiciously elevated recombination
#     rates — typically sub-ceiling bad markers that remained after round 1.
#     The v2 file is the output of step 12 (HMM re-run with round-1 markers
#     masked), so the round-1 bad markers are already gone.
#
#   This script auto-detects which round to run based on whether a v2 file
#   is present. If v2 is found it runs round 2; otherwise it runs round 1
#   on the v1 file.
#
#   This is a BACKCROSS (BC) population. Expected genotype ratios depend on
#   which parent is the recurrent parent:
#     AA-recurrent BC:  AA~50%,  AB~50%,  BB~0%
#     BB-recurrent BC:  AA~0%,   AB~50%,  BB~50%
#
#   The script auto-detects BC direction by sampling genotype frequencies at a
#   known-good marker (r < 0.05, well below any threshold). If AA >> BB it is
#   AA-recurrent; if BB >> AA it is BB-recurrent.
#
#   Each flagged marker is scored:
#     score = NN% + |AA% - exp_AA| + |AB% - 50| + |BB% - exp_BB|
#   The marker with the higher score is the bad one.
#   If a marker is sandwiched between two consecutive flagged intervals
#   (shared marker rule), it is flagged regardless of score, with a WARNING
#   printed if the badness score disagrees.
#
# USAGE:
#   bash check_bad_markers.sh
#   (no arguments — files are auto-detected from the current directory)
#
# OUTPUT:
#   - Prints genotype frequencies and badness scores for each flagged interval
#   - MERGES newly identified bad markers into bad.marks.txt (does not
#     overwrite — existing entries from round 1 are preserved)
#   - Format: chromosome <tab> position, no header, sorted
# =============================================================================

# Recombination rate thresholds:
#   Round 1 (v1 file): 0.25 — catches markers where the HMM hit the MaxRR
#                              ceiling, a hard sign of a bad marker
#   Round 2 (v2 file): 0.1  — catches sub-ceiling bad markers that survived
#                              round 1 after the obvious bad markers are masked
RR_THRESHOLD_ROUND1=0.2
RR_THRESHOLD_ROUND2=0.1

BADMARKS="bad.marks.txt"

# Auto-detect intrascaff file: prefer v2 (round 2) over v1 (round 1).
# The threshold is set automatically based on which file is found.
INTRASCAFF=$(ls *intrascaff*v2*.txt 2>/dev/null | head -1)
if [[ -n "$INTRASCAFF" ]]; then
    RR_THRESHOLD=$RR_THRESHOLD_ROUND2
    echo "Found v2 intrascaff file — running ROUND 2 (threshold = $RR_THRESHOLD)"
else
    INTRASCAFF=$(ls *intrascaff*v1.txt 2>/dev/null | head -1)
    RR_THRESHOLD=$RR_THRESHOLD_ROUND1
    echo "No v2 file found — running ROUND 1 (threshold = $RR_THRESHOLD)"
fi
if [[ -z "$INTRASCAFF" ]]; then
    echo "ERROR: No intrascaff file found matching *intrascaff*v2*.txt or *intrascaff*v1.txt"
    exit 1
fi
echo "Using intrascaff file: $INTRASCAFF"

# Auto-detect the individual list — looks for any file matching *indivi_list.txt
INDIVI_LIST=$(ls *indivi_list.txt 2>/dev/null | head -1)
if [[ -z "$INDIVI_LIST" ]]; then
    echo "ERROR: No individual list found matching *indivi_list.txt in current directory"
    exit 1
fi
echo "Using individual list: $INDIVI_LIST"
echo ""

# Count total individuals — used later to compute percentages
TOTAL_INDIVI=$(wc -l < "$INDIVI_LIST")
echo "Total individuals: $TOTAL_INDIVI"
echo ""

# -----------------------------------------------------------------------------
# STEP 1: Auto-detect BC direction (AA-recurrent vs BB-recurrent)
#
# In a backcross, one homozygous class should be ~50% and the other ~0%.
# We sample the first non-flagged marker in the intrascaff file (r < 0.2499)
# and count AA vs BB calls across all individuals.
# If AA > BB by a clear margin → AA-recurrent → expected: AA=50, AB=50, BB=0
# If BB > AA by a clear margin → BB-recurrent → expected: AA=0,  AB=50, BB=50
# -----------------------------------------------------------------------------

# Find the first marker with a recombination rate clearly below any threshold
# we might use (r > 0 and r < 0.05). Using 0.05 as the cutoff here ensures
# the reference marker is unambiguously good regardless of RR_THRESHOLD.
ref_scaff=""
ref_pos=""
while IFS=$'\t' read -r scaffold pos rate; do
    is_good=$(awk -v r="$rate" 'BEGIN { print (r+0 > 0 && r+0 < 0.05) ? 1 : 0 }')
    if [[ $is_good -eq 1 ]]; then
        ref_scaff=$scaffold
        ref_pos=$pos
        break
    fi
done < "$INTRASCAFF"

# Count AA and BB at the reference marker across all individuals
ref_aa=0; ref_bb=0
while read -r plantID; do
    gfile="g.Genotypes.${plantID}.bam.txt"
    [[ ! -f "$gfile" ]] && continue
    call=$(awk -F'\t' -v sc="$ref_scaff" -v pos="$ref_pos" \
        '$2==sc && $3==pos {print $4; exit}' "$gfile")
    case "$call" in
        AA) ((ref_aa++)) ;;
        BB) ((ref_bb++)) ;;
    esac
done < "$INDIVI_LIST"

# Determine BC direction: if AA is more than 3x BB → AA-recurrent (and vice versa)
bc_dir=$(awk -v aa=$ref_aa -v bb=$ref_bb 'BEGIN {
    if (aa > 3*bb) print "AA"
    else if (bb > 3*aa) print "BB"
    else print "UNKNOWN"
}')

if [[ "$bc_dir" == "AA" ]]; then
    exp_AA=50; exp_AB=50; exp_BB=0
    echo "BC direction detected: AA-recurrent (AA=$ref_aa, BB=$ref_bb at $ref_scaff:$ref_pos)"
    echo "Expected ratios: AA~50%  AB~50%  BB~0%"
elif [[ "$bc_dir" == "BB" ]]; then
    exp_AA=0; exp_AB=50; exp_BB=50
    echo "BC direction detected: BB-recurrent (AA=$ref_aa, BB=$ref_bb at $ref_scaff:$ref_pos)"
    echo "Expected ratios: AA~0%  AB~50%  BB~50%"
else
    # Cannot determine direction — fall back to F2-like midpoint as a neutral score
    exp_AA=25; exp_AB=50; exp_BB=25
    echo "WARNING: BC direction unclear (AA=$ref_aa, BB=$ref_bb at $ref_scaff:$ref_pos)"
    echo "Falling back to symmetric expected ratios: AA~25%  AB~50%  BB~25%"
    echo "Consider checking your population type."
fi
echo ""

# -----------------------------------------------------------------------------
# STEP 2: Find all intervals where r = 0.25 (hit the MaxRR ceiling)
# We read the intrascaff file line by line, keeping track of the previous line
# so we always know the pair: marker_above (prev line) and marker_below (curr line)
# The recomb rate on a line is the rate BETWEEN that marker and the NEXT marker,
# so when we see r=0.25 on line j, the bad interval is between marker j and j+1
# -----------------------------------------------------------------------------

echo "=========================================="
echo " Flagged intervals with r >= $RR_THRESHOLD"
echo "=========================================="

# We'll collect the flagged marker pairs into arrays for genotype checking below
declare -a SCAFFOLDS
declare -a MARKER_ABOVE   # marker j   (has r=0.25 in output)
declare -a MARKER_BELOW   # marker j+1 (next line after r=0.25)

prev_scaffold=""
prev_pos=""
prev_rate=""
next_line_is_below=0
pair_idx=0

while IFS=$'\t' read -r scaffold pos rate; do

    # If the previous line had r=0.25, this line is the "marker below" the bad interval
    if [[ $next_line_is_below -eq 1 ]]; then
        SCAFFOLDS[$pair_idx]=$prev_scaffold
        MARKER_ABOVE[$pair_idx]=$prev_pos
        MARKER_BELOW[$pair_idx]=$pos
        echo "  Interval $((pair_idx+1)): $prev_scaffold  above=$prev_pos  below=$pos  (r=$prev_rate)"
        ((pair_idx++))
        next_line_is_below=0
    fi

    # Check if this line's rate exceeds RR_THRESHOLD
    # Use awk for floating point comparison (bash can't do floats)
    is_bad=$(awk -v r="$rate" -v thr="$RR_THRESHOLD" 'BEGIN { print (r+0 >= thr+0) ? 1 : 0 }')
    if [[ $is_bad -eq 1 ]]; then
        next_line_is_below=1
    fi

    prev_scaffold=$scaffold
    prev_pos=$pos
    prev_rate=$rate

done < "$INTRASCAFF"

echo ""

# -----------------------------------------------------------------------------
# STEP 3: For each flagged interval, count genotype frequencies at both markers
# and compute a "badness" score for each.
#
# Score = NN% + |AA% - exp_AA| + |AB% - 50| + |BB% - exp_BB|
# A perfect BC marker scores 0. A bad marker scores high.
# Expected values are set based on the auto-detected BC direction above.
#
# We also track which markers appear in more than one interval — a marker
# shared between two consecutive 0.25 intervals is sandwiched by bad signal
# on both sides. This is used as a confirmation check against the badness score.
# -----------------------------------------------------------------------------

declare -a SCORE_ABOVE
declare -a SCORE_BELOW

for (( i=0; i<pair_idx; i++ )); do

    scaff=${SCAFFOLDS[$i]}
    above=${MARKER_ABOVE[$i]}
    below=${MARKER_BELOW[$i]}

    echo "=========================================="
    echo " Interval $((i+1)): $scaff"
    echo " Marker ABOVE interval: position $above"
    echo " Marker BELOW interval: position $below"
    echo "=========================================="

    # For each of the two markers, count genotypes across all individuals
    for marker_pos in "$above" "$below"; do

        aa=0; ab=0; bb=0; nn=0

        while read -r plantID; do
            gfile="g.Genotypes.${plantID}.bam.txt"
            if [[ ! -f "$gfile" ]]; then
                continue
            fi

            # Extract the genotype call for this scaffold + position from the individual's file
            # Column layout: plantID  scaffold  position  genotype
            call=$(awk -F'\t' -v sc="$scaff" -v pos="$marker_pos" \
                '$2==sc && $3==pos {print $4; exit}' "$gfile")

            case "$call" in
                AA) ((aa++)) ;;
                AB) ((ab++)) ;;
                BB) ((bb++)) ;;
                NN) ((nn++)) ;;
            esac

        done < "$INDIVI_LIST"

        # Calculate percentages — tells us whether this marker looks like a
        # healthy BC marker or something aberrant
        total_called=$((aa + ab + bb + nn))
        pct_aa=$(awk -v n=$aa -v t=$total_called 'BEGIN { printf "%.1f", (t>0) ? n/t*100 : 0 }')
        pct_ab=$(awk -v n=$ab -v t=$total_called 'BEGIN { printf "%.1f", (t>0) ? n/t*100 : 0 }')
        pct_bb=$(awk -v n=$bb -v t=$total_called 'BEGIN { printf "%.1f", (t>0) ? n/t*100 : 0 }')
        pct_nn=$(awk -v n=$nn -v t=$total_called 'BEGIN { printf "%.1f", (t>0) ? n/t*100 : 0 }')

        # Compute badness score using BC-appropriate expected ratios
        # Higher score = greater deviation from what a healthy BC marker should look like
        score=$(awk -v aa=$pct_aa -v ab=$pct_ab -v bb=$pct_bb -v nn=$pct_nn \
                    -v eaa=$exp_AA -v eab=$exp_AB -v ebb=$exp_BB \
            'BEGIN {
                score = nn + sqrt((aa-eaa)^2) + sqrt((ab-eab)^2) + sqrt((bb-ebb)^2)
                printf "%.1f", score
            }')

        echo "  Position $marker_pos:"
        echo "    AA: $aa  ($pct_aa%)   AB: $ab  ($pct_ab%)   BB: $bb  ($pct_bb%)   NN: $nn  ($pct_nn%)   badness score: $score"

        # Store scores for later comparison
        if [[ "$marker_pos" == "$above" ]]; then
            SCORE_ABOVE[$i]=$score
        else
            SCORE_BELOW[$i]=$score
        fi

    done

    echo ""
    echo "  --> Expected healthy BC marker: AA~${exp_AA}%  AB~${exp_AB}%  BB~${exp_BB}%  NN~low"
    echo "  --> The marker deviating most from these ratios is likely the bad one"
    echo ""

done

# -----------------------------------------------------------------------------
# STEP 4: Determine which marker to remove for each flagged interval
#
# Rule 1 — Isolated single interval: badness score decides.
#   The marker with the higher score is removed.
#
# Rule 2 — Back-to-back intervals: shared marker rule decides.
#   The marker sandwiched between two consecutive 0.25 intervals is removed.
#   Badness scores are still computed for both flanking markers and used as a
#   confirmation check. If the shared marker has a BETTER badness score than
#   the outer marker, a WARNING is printed to stdout (visible in the Slurm
#   .out file) so you can manually inspect and correct if needed.
#
# Deduplicate the final list — if the same marker is flagged by multiple
# intervals, write it to bad.marks.txt only once.
# -----------------------------------------------------------------------------

echo "=========================================="
echo " Bad marker decisions"
echo "=========================================="

# Use an associative array to deduplicate: key = "scaffold<tab>position"
declare -A BAD_MARKERS

for (( i=0; i<pair_idx; i++ )); do

    scaff=${SCAFFOLDS[$i]}
    above=${MARKER_ABOVE[$i]}
    below=${MARKER_BELOW[$i]}
    score_a=${SCORE_ABOVE[$i]}
    score_b=${SCORE_BELOW[$i]}

    # Check if this marker is shared with an adjacent interval:
    # - "above" of interval i == "below" of interval i-1 → above is the shared bad marker
    # - "below" of interval i == "above" of interval i+1 → below is the shared bad marker
    shared_above=0
    shared_below=0

    if (( i > 0 )) && [[ "${MARKER_BELOW[$((i-1))]}" == "$above" && "${SCAFFOLDS[$((i-1))]}" == "$scaff" ]]; then
        shared_above=1
    fi
    if (( i < pair_idx-1 )) && [[ "${MARKER_ABOVE[$((i+1))]}" == "$below" && "${SCAFFOLDS[$((i+1))]}" == "$scaff" ]]; then
        shared_below=1
    fi

    if [[ $shared_above -eq 1 || $shared_below -eq 1 ]]; then

        # --- Back-to-back: shared marker rule decides ---
        if [[ $shared_above -eq 1 ]]; then
            bad_marker="$above"
            outer_score=$score_b
            bad_score=$score_a
        else
            bad_marker="$below"
            outer_score=$score_a
            bad_score=$score_b
        fi

        # Confirm with badness score — warn if they disagree
        scores_agree=$(awk -v bs="$bad_score" -v os="$outer_score" \
            'BEGIN { print (bs+0 >= os+0) ? 1 : 0 }')

        if [[ $scores_agree -eq 1 ]]; then
            reason="shared marker rule (confirmed by badness score: $bad_score vs $outer_score)"
        else
            reason="shared marker rule (badness score: $bad_score vs $outer_score)"
            echo ""
            echo "  *** WARNING: interval $((i+1)) $scaff — shared marker rule removes $bad_marker"
            echo "  *** but its badness score ($bad_score) is BETTER than the outer marker ($outer_score)"
            echo "  *** Please inspect this interval manually before accepting this decision"
            echo ""
        fi

    else

        # --- Isolated interval: badness score decides ---
        is_above_worse=$(awk -v sa="$score_a" -v sb="$score_b" \
            'BEGIN { print (sa+0 >= sb+0) ? 1 : 0 }')
        if [[ $is_above_worse -eq 1 ]]; then
            bad_marker="$above"
            reason="higher badness score ($score_a vs $score_b)"
        else
            bad_marker="$below"
            reason="higher badness score ($score_b vs $score_a)"
        fi

    fi

    echo "  Interval $((i+1)): $scaff — remove position $bad_marker ($reason)"

    # Add to deduplicated bad markers list (tab-delimited key)
    BAD_MARKERS["${scaff}"$'\t'"${bad_marker}"]=1

done

echo ""

# -----------------------------------------------------------------------------
# STEP 5: Merge newly identified bad markers into bad.marks.txt
#
# Format: chromosome <tab> position  (no header, tab-delimited)
# This file is read directly by step 12 to mask bad markers before re-running
# the HMM.
#
# MERGE behaviour: existing entries in bad.marks.txt (from round 1) are
# preserved. New markers identified here are added. Duplicates are removed.
# This ensures bad.marks.txt always contains the complete list across all
# rounds.
# -----------------------------------------------------------------------------

# Load existing bad markers into the associative array so they are preserved
if [[ -f "$BADMARKS" ]]; then
    while IFS=$'\t' read -r scaff pos; do
        [[ -z "$scaff" ]] && continue
        BAD_MARKERS["${scaff}"$'\t'"${pos}"]=1
    done < "$BADMARKS"
    echo "Existing bad markers loaded from $BADMARKS: ${#BAD_MARKERS[@]} entries"
fi

# Write complete merged list
> "$BADMARKS"
for key in "${!BAD_MARKERS[@]}"; do
    echo -e "$key" >> "$BADMARKS"
done

# Sort by chromosome then position for readability
sort -k1,1 -k2,2n "$BADMARKS" -o "$BADMARKS"

echo "=========================================="
echo " Final bad.marks.txt (all rounds merged):"
cat "$BADMARKS"
echo "=========================================="
