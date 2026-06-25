# Segregation_Distortion_Mapping
Follow: dx.doi.org/10.17504/protocols.io.bjnbkman for processing ddRAD data from raw fastqs to vcf (2021). The pipeline we have here is modified from it.
Make sure to request **i5** and **i7** fastq files from the sequencing center.
## Step 1: remove PCR duplicates from each library
Run: `step1_rmdup.sh`

This step removes PCR duplicates using the i5 molecular barcode.

Python script: `step1_rmdup.py` (requires unzipped fastq files), runs in about 2 hours per library.

**Required inputs:**
R1, R2, i5 and i7 fastq files
Rename files from the sequencing center to match the format of R1, R2, i5, and i7 fastq files.
Files from YCGA will be supplied as:
- I1 = I7 index
- R1 = Sequencing Read 1
- R2 = I5 index read
- R3 = Sequencing Read 2

**⚠️Note:** This Python script will not recognize your fasta files unless their name is ${library}_R1.fastaq. The prefix MUST be the library name, and the suffix MUST be R1, R2, i7, or i5. Otherwise, it won't input it.

**To prefrom this step:** run `step1_rmdup.sh` script. This bash script will call the Python script `step1_rmdup.py`, so make sure the Python script is in the same directory as the bash file.
This step will remove duplicates and output fastq.gz files for i5, i7, R1, and R2
The resulting files will have the same prefix, but will have the suffix .rmdup.1.fastq (forward reads) and .rmdup.2.fastq (reverse reads).

(skipping the step to flip the reads from the Fishman lab protocol, and moving straight to the next step to demultiplex samples in each library)

## Step 2: Sample Demultiplexing using Stacks
Run:  `step2_demultiplex.sh`

Requires a .txt file with sample barcodes of each well on the plate to demultiplex samples.
Prepare a txt file listing sample IDs and corresponding barcodes. The barcodes can be found in `A1_A1_barcodes`, note that the sample IDs will be different for rach library.

**⚠️Important Note:** The BestRAD protocol we used to construct the library generates a unique “GG” at the beginning, so you have to add “GG” before you formal barcodes.
 
**⚠️Important Note:** Before running this script, rename fastq.gz files to a name format as received from the sequencing center. Otherwise, stacks will not recognize the input files
rename to format PS_P4_L1_2024_S3_L006.R1_001.fastq.gz (same for R2 file) or PS_P4_L1_2024_S3_L006_R1_001.fastq.gz (or with underscore before R#)

Stacks is picky about the length of the file names for each sample. picky how many "columns" it can have (meaning how many sets of characters separated by _ (did not test if the issue is the overall number of characters or if _ creates a new "column" to the file)
The barcodes file needs to be a tab-delimited .txt file.

## Step 3: Adaptor Trimming
Run: `step3_trim.sh`

software: trimmomatic

input:
4 fastq files per sample generated from Stacks
`step3_adapter.fa` file contain adapter sequences to be trimmed

Details for the adapter sequences file (generated with a text editor):
- Adapter 1 & 2: forward and reverse adapter sequences from the BestRAD oligos (including the GG before the formal barcodes).
- Adaptor 3 & 4: forward and reverse adapter sequences from NEBNext adapters kit as listed in [NEBNext Primer instruction manual](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7335_e7500_-e7710_e7730.pdf?rev=2e735fd18b544d46b36ee0e88353ef5c&sc_lang=en-us&hash=CC77B45817715F3ED3A8F3B1953450EB)
- Forward and reverse i5 adapter sequences.
- Forward and reverse i7 adapter sequences (modified based on i7 adapter sequence, nucleotide number, might vary between 6-8).

Note that the NNNNNN in the sequences is the unique identifier for each sample/individual.
## Step 4: Sequence alignment, Sorting, Cleaning, and Indexing
Run: `step4_align.sh` 

Software: `BWA` and `SAMtools`

Input: `fastq.gz` files, the output from trimmomatic, and reference genome (indexed with samtools faidx). Then the SAMtools uses the sam files generated from BWA to sort and clean the sequences.

Output: `.bam files`

Sorted files are smaller in size and faster to process; this is why the downstream tools require sorted files. The cleaning step removes aligned sequences with low scores.

Our lab already has the index for the reference genome in the same directory, so don't worry about making an index.
## Step 5: SNP Calling
Run: `step5_VCF.sh` and `bam_list_loop.sh` to make the bam list.

Software: `BCFtools` and `SAMtools`

Input: a .txt list of full paths to all the bam files needed for this VCF. I have created a loop that can make this list for you. Also, make sure to add the parents and F1. Each VCF should be a linkage group (do not mix individuals from different cross directions unless there is a reason to do so, which I might do in the future)

**⚠️Note:** This script DOES NOT filter the VCF and keeps multiallelic variants (i.e., does not filter to only keep biallelic sites). This is **important** as VCF filtering will take place downstream when calling ancestry is step 6. I also added two extra lines to create an index and a summary statistics file. The reason why this script does not do any filtering is that it can create a "raw" VCF that then can be filtered to different software depending on their requirements. 

## Step 6: VCF processing using 3 Python scripts
Run: `step6_VCF_prog1.py`, `step6_VCF_prog2.py`, and `step6_VCF_prog3.py` back to back.
- **step6_VCF_prog1.py:** This script processes a compressed VCF file containing SNP calls and filters variants based on quality and allele balance criteria. It reads a .vcf.gz file for a specified library (e.g., "1A"), extracts SNPs with a minimum mapping quality (MQ ≥ 20), and evaluates allele depth (AD) across samples. For each SNP, it calculates the number of samples with valid calls, average read depth, and average reference allele proportion. SNPs that pass thresholds for sample coverage (≥ 50 samples) and allele balance (between 0.2 and 0.8) are written to a slimmed-down VCF and two read depth summary files. It also tracks the number of SNPs per scaffold and the position of the last SNP, outputting this to a scaffold summary file. You can change the: Min_MQ_score, Min_lines_called, minQ (allele freq),maxQ(allele freq).

**⚠️Note:** for my own crossing design (BC and not F2s), I have changed the filtering ratios since they will change from one direction to another based on who is the recurrent parent in a given backcrossing population. See the example below:

```python

#-------------------------------------------------------------------------------
# Step6 Prog1 — VCF filtering for map 7C and 8C
# BC to CCC9 = BB parent
# qR near 0.0 = strong BB excess (distortion); qR near 0.5 = all heterozygotes
# upper bound 0.70 excludes AA-biased sites impossible in a BC to BB
#   minQ = 0.05
#   maxQ = 0.70
# Change the script to the values below when the recurrent parent is not IM62 or what the genome is aligned too

  Min_MQ_score = 20
    Min_lines_called = 50
    minQ = 0.05
    maxQ = 0.70

#-------------------------------------------------------------------------------
# Step6 Prog1 — VCF filtering for map 5C and 6C
# BC to IM62 = AA parent
# qR near 0.5 = all heterozygotes; qR near 1.0 = strong AA excess (distortion)
# lower bound 0.30 excludes BB-biased sites impossible in a BC to AA
#   minQ = 0.30
#   maxQ = 0.95
#-------------------------------------------------------------------------------
    Min_MQ_score = 20
    Min_lines_called = 50
    minQ = 0.30
    maxQ = 0.95
```

-  **step6_VCF_prog2.py:** This is a diagnostic script. Run this script two times. First run produced info about the depth to decide on the cut-offs. Reads with more than 15 are usually at repetitive genomic regions. Lower the cut-off to produce the final SNP.limited.txt file when you run the script again. It is okay not to be super strict here since we will do depth filtering in step9B at a window level. I did the first run at 50, then decided that 15 was a good number to remove most outliers.
**⚠️Important Note:** Make sure to look at the .out file from the first Python program and see which line contains the two parents and use this info to edit the second Python script (ln 62). Also, don't forget to change the number of individuals.

-  **step6_VCF_prog3.py:** This script processes a list of selected SNPs (SNPs.limited.txt) and extracts reference and alternate allele depths for each individual across those SNPs, organizing the output into scaffold-specific and sample-specific files. It will output a .txt file for all .bam files in their own directories, along with a .txt file for each chromosome. Also, it should output an `all.samples.txt` file. I have changed this script to output the files to a new directory called "refAlt". 

**⚠️Important Note:** that the third Python script outputs lots of data (.bam.txt file for each bam file in the original aligned directory) AND does not overwrite the generated files if it fails and is re-run. Instead, it will just keep adding lines to the pre-existing files. Thus, always make sure to remove any generated files if your script failed to run before re-running again.

## Step 7: Assign parental ancestry (R)
Run: `step7_assign_parental_ancestry.R` as a cluster job, as it will take a few hours to run. Use `step7_run_R_parental_ancestry.sh` to run. 

This R script compares allele depths from two parental lines (e.g. IM62 and IMPO) to classify each SNP site as REF (IMPO), ALT(IM62), HET, or missing (NN), then uses this classification to polarize genotypes in a set of recombinant individuals, flipping allele counts where necessary to align with parental ancestry. For the polarization step, if parents are opposite (e.g., IMPO = REF, IM62 = ALT), it retains allele counts. - If parents are flipped (IMPO = ALT, IM62 = REF), swaps ref and alt counts to align with ancestry.

I have modified the R script to exclude the two parental files so I don't corrupt them in case this step fails. Second, after this script and onward, AA = CCC9 or IMPO and BB = IM62

Output: `parental_ancestry.tt`, update the `.nam.txt' file to make sure they are poralized. and `all.scaffold_*txt.`

## step 8: Windows
Run: `step8_window.py`. Use `step8_window_run.sh`

This step averages across a window based on the number of SNPs/ number of reads/ some range. You can define heterozygote calls as het deviation (I changed it to 0.2 to make it stricter when calling heterozygous windows) and the minimum numebr of SNPs per window. I choose to relax the SNPs/window filter here as I will filter in step 9 for that. This will output 3 files for each individual. 1) a Genotypes file, 2) a genostats file (total number of windows with AA, AB, NN genotypes), and 3) a windows file (which I think is a count for each site that was used to build the genotypes file).
**⚠️Important Note:** I hated how this script outputs files with no headers, so I added a few Python lines to add headers to the output files. However, downstream, these headers can cause an issue, which I tried to fix in the downstream script. For now, I added a # in front of those lines to just avoid the headache altogether. 

**⚠️Important Note:** Future Hagar noticed that when there are markers with high percentage of NN, this inflates the peak and cause an artifcat downstream. one way to delet these windows is to make sure to add them possibly (site_list) to be exlcuded or include them into the bad.marks.txt file? 

The goal: identify genomic windows where more than 10% of individuals have genotype = NN (missing).

```bash
# Run from the directory containing Genotypes.*.txt for a given map
cat Genotypes.*.txt | \
  awk 'NR==1 || $4=="Genotype" {next}          # skip headers
       {key=$2"\t"$3; total[key]++; 
        if($4=="NN") nn[key]++}
  END {for(k in total)
         print k"\t"(nn[k]+0)"\t"total[k]"\t"(nn[k]+0)/total[k]}' | \
  awk '$4 > 0.10 {print $1"\t"$2}' \
  > badmarks_highNN.txt
```

This produces a Scaffold\tWindowStart file — exactly the format step12 expects. Adjust the 0.10 threshold as needed (10% is conservative; you might use 0.15–0.20 depending on how strict you want to be).

There is a filtering step at step 9B, but I will investigate further. 


## Step 9: Second Ancestry loop to filter and output .g files (R)
This step will process the genotypes, genostats, and windows output from step 8. It will filter the genotype files and also filter out windows with low depth.

💡 I split this step into 9A and 9B. 9A takes a little bit of time to generate the site counts, it also detect the window number automatically for each map. 9B filters bad individuals and windows-based and generates the g. files. Check and change the parameters at the top of the script. 


Excluded any individuals that:

**(1)** have 80% or of their windows missing or (NN)

**(2)** are 70% homozygous for the backcrossed parent or 70% heterozygous.

**(3)** are 10% homozygous for the non-recurrent parent (this can be tighter).


Exclude any windows that:

**(1)** have less than 3 SNPs to avoid allele dropout.

**(2)** absent in more than 90% of the population (i.e. very few individuals carry that marker).

Step 9B will output the list of bad individuals and sites then only generate g. files for the passes individuals with passed windows. 


**⚠️Note:**  The Python script from step 8 will output files with a header if you decided to use the header line in the script, so make sure that header = TRUE if that was the case.
 

## Step 10 GOOGA: Galculate Genotype Error Rate

Run: `step10_calc_genotype_err_rate.py` as a parallel job using this job list `step10_genotypes_list.txt`. This runs for each sample; use `step10_sample_list_loop.sh` to make the list.

**⚠️Important Note:** make sure to make an empty file called `bad.marks.txt` or the script won't run
**⚠️Important Note:** Make sure the genotype ratios are correct based on the cross direction (e.g. start_probability = {'AA':0.1,'AB':0.4,'BB':0.5}). Notice that I allowed for a little bit of AA windows just so I don't get errors. No need to change the transition probabilities for each cross direction; I have changed this to work with a BC design, and it's symmetrical, so no need to change further. 
```bash
module load dSQ
dsq --job-file step10_genotypes_list.txt --mem-per-cpu 4g -t 20:00 --mail-type ALL
#then run the sbatch code it will generate
```
Note: this Python script cannot run with barebone Python 2, this is why I am loading `SciPy-bundle/2020.11-foss-2020b-Python-2.7.18`

This Python script will generate a file for each .g file, and each file will contain one line. I used a loop called `combine_genotype_err_loop.sh` to combine them in a txt file called `1A_combined_genotypes_err.txt`. In this step, you will notice that there are errors in calculating the genotype error rate for some individuals. This happens, I believe, due to big gaps in some of the chromosomes that give that error. I have modified that Python script to prevent this error. But I still had some problematic individuals (just one or two) I exlcuded them and made a new list named `1A_filtered_genotypes_err_rates.txt`.

In Excel, filter individuals based on their genotype errors. John's criteria is less than 20% (.2) error for each of the different error rates. The columns should be: Sample ID, Marker number, e1, e2, beta, and likelihood. In Excel, make a new columns called drop? nad se the AND() for e1, e2, and beta `=AND(C2<0.2,D2<0.2,E2<0.2)` tell you which samples are high quality - i.e. any individuals that have <.2 for all error rates. Save as tab-delimited file,

💡 To automate doing this with a script, run `combine_and_filter_genotype_err.sh` and it will generate a filtered list of genotype error and a list containing the names of these individuals used in step 11.


## Step 11 GOOGA: Calculate intra-scaffold recombinational fractions
Run: `step11_hmm.intrascaff.R.py` using `step11_intrascaff_rec_rates.sh`. This will take a few hours.

Input: `1A_filtered_genotypes_err_rates.txt`(has to be tab-delimited), a random `g.Genotypes.PlantID.bam.txt` to get the markers names, `1A_low_genotype_err_list.txt`, which is just a list of the individuals in the genotype error rates file (make sure it's the naked name with no `.bam.txt`. Finally, an empty `bad.marks.txt` 

 **⚠️Important Note:** In `1A_filtered_genotypes_err_rates.txt`, the columns must be tab-delimited or else the script won't run.

 **⚠️Important Note:** Make sure that in the Python script `step11_hmm.intrascaff.R.py` the AA, AB, and BB ratios and transition probabilities are reflecting the correct ratios depending on the parent F1s that were backcrossed to. This will also depend on your definetion on which parent is AA and which is BB in step 7
 
The output file (1A_intrascaff_v1.txt) has your first linkage map! Column 1 is the chromosome, column 2 is the marker. In column 3, each value is the intrascaffold recombination rate between the marker for that row and the
marker below. The units on this are in probability of recombination (or cM/100). The maximum value here is 0.25 (this should not be a surprise because the max recombination rate is .5 (ie. unlinked) and there is a 50% chance that recombination occurs before the marker and a 50% chance it occurs after the marker (thus .25 each way). This may mean that a marker is either in the wrong location(in need of GOOGA!) or just poorly fits the genotype error algorithm and should have been discarded.

How do you determine which of these possibilities is correct? Well, we can discard these markers and see if there are still large recombination distances in the same spots. If it is the case, then our
chromosome needs reordering in Googa. If this fixes the problem, then the marker was simply bad.

Ok, let's get rid of the bad markers. So, highlight each of the crappy markers in yellow in Excel. There is likely either one maker with .25 or two makers. Remember, the recombination rate is calculated with the marker below it. So get rid of the marker below the .25 in case of a single marker with a .25. If there are two markers back to back with .25, then get rid of the second one. Compile the list of bad markers into your previously created tab-delimited text file named `bad.marks.txt`. You only need two columns, Chromosome and Marker, but do not put in a header line. Now we can rerun the program calculating intrascaffold error rates with the `bad.marks` version of the last program and see if we get rid of the 0.25 errors. It is normal for the last marker in the chromosome to have a large negative value, so don't delete that. 


💡 I created a script that automatically checks bad markers and outputs them in the `bad.marks.txt`. run `check_bad_marks.sh` script. This script can be run twice; the first run auto-detects `intrascaff_v1.txt`. Discover the bad markers, then when step 12 is complete, re-run `check_bad_marks.txt` and it will auto detect `intrascaff_v2.txt` and apply a stricter threshold (i.e. going from 0.25 to 0.1). This script is also better than doing the bad marker removal manually, as it traces back each window and calculates a score based on which it chooses which window to exclude and which to ket (i.e. not always the second window)

## Step 12: GOOGA: Calculate Rntra-scaffold Recombinational Fractions After Removing Bad Markers
Run: `step12_hmm.intrascaff.R.badmarks.py` using `step12_intrascaff_rec_badmarks.sh`. This will take a few hours.

Input: `1A_filtered_genotypes_err_rates.txt`, a random genotype `g.Genotypes.1A_1A_D8.bam.txt`, a list of individuals "1A_low_genotype_err_list.txt`, and `bad.marks.txt`.

Output: `1A_intrascaff_v2.txt`

## Step 13: Estimate Genotypes at Each Locus for Each Individual
Before running this Python script, we will need to calculate the cumulative recombination fraction for each marker for each scaffold/chromosome. I used Excel to do so using the following instructions:
- First marker, cM = 0
- Second marker, cM = 100*D2
- Third maker, cM= 100*SUM($D$2:D3), then drag till the end of the chromosome
- Assuming your recombination fractions are in column D, with the first marker starting in cell D2. Start each chromosome over with this same code and then fill down to the end of the chromosome with the formula for the 3rd marker. Once done, copy and paste values only into a new sheet to avoid missing values. The final file should look like `1A_intrascaff_v3.txt`, it must be tab-delimited. The wrong format of this file can cause errors.

💡 I have automated this step; run `make_intrascaff_v3.sh`

Run: `step13_run_genotype.pp.txt` that uses `step13_genotype.pp.py`. This should run as a parallel job on the cluster for each linkage group (i.e. Chr_01..14)

```bash
module load dSQ
dsq --job-file step13_run_genotype.pp.txt --mem-per-cpu 4g -t 30:00 --mail-type ALL
#then run the sbatch code it will generate
```

Output: This script outputs 14 files for each chromosome with posterior probabilities for each marker. `Chr_01..Chr_14.pp.txt`

I then use the following code to concatenate all of them. Make sure that all chromosomes are present after concatenating, because the length of this file (i.e. the number of rows) should equal the number of markers in `1A_intrascaff_v3.txt`
```bash
cat Chr_{01..14}.pp.txt > 1A_allmarkers.txt
```

Then we move to R to make our genotype matrix and use rQTL to plot our data!!!!

The R script is called `1A_googa_processing.Rmd`. I usually do this on the cluster's R



