# Segregation_Distortion_Mapping
follow: dx.doi.org/10.17504/protocols.io.bjnbkman for processing ddRAD data from raw fastqs to vcf (2021).
Make sure to request **i5** and **i7** fastq files from the sequencing center.
## Step 1: remove PCR duplicated from each library
Run: `step1_rmdup.sh`

This step removes PCR duplicate using the i5 molecular barcode.

python script: `step1_rmdup.py` (requires unzipped fastq files), runs in about 2 hours per library.

**Required inputs:**
R1, R2, i5 and i7 fastq files
Rename files from sequencing center to match format of R1, R2, i5 and i7 fastq files:
files from YCGA will be supplied as:
- I1 = I7 index
- R1 = Sequencing Read 1
- R2 = I5 index read
- R3 = Sequencing Read 2

**To prefrom this step:** run `step1_rmdup.sh` script. This bash script will call the python script `step1_rmdup.py`.
This step will remove duplicates and output fastq.gz files for i5, i7, R1 and R2
The resulting files will have the same prefix, but will have the suffix .rmdup.1.fastq (forward reads) and .rmdup.2.fastq (reverse reads).

(skipping the step to flip the reads from the Fishman lab protocol, moving straight to next step to demultiplex samples in each library)

## Step 2: Sample Demultiplexing using Stacks
Run:  `step2_demultiplex.sh`

Requires a .txt file with sample barcodes of each well on the plate to demultiplex samples.
Prepare a txt file listing sample IDs and corresponding barcodes. The barcodes can be found in `A1_A1_barcodes`, note that the sample IDs will be different for rach library.

Important: The BestRAD protocol we used to construct the library generates a unique “GG” at the beginning, so you have to add “GG” before you formal barcodes.

**⚠️Important** before running this script, rename fastq.gz files to name format as received from the sequencing center. Otherwise stacks will not recognize the input files
rename to format PS_P4_L1_2024_S3_L006.R1_001.fastq.gz (same for R2 file) or PS_P4_L1_2024_S3_L006_R1_001.fastq.gz (or with underscore before R#)

stacks is picky about the length of the file names for each sample. picky how many "columns" it can have (meeaning how many sets of characters separated by _ (did not test if the issue is the overall number of characters or if _ creates a new "column" to the file)
The barcodes file needs to be a tab delimited .txt file.

## Step 3: Adaptor Trimmning
Run: `step3_trim.sh`

software: trimmomatic

input:
4 fastq files per sample generated from stacks
`step3_adapter.fa` file contains adapter sequences to be trimmed

Details for the adapter sequences file (generate with text editor):
- Adapter 1 & 2: forward and reverse adapter sequences from the BestRAD oligos (including the GG before the formal barcodes).
- Adaptor 3 & 4: forward and reverse adapter sequences from NEBNext adapters kit as listed in [NEBNext Primer instruction manual](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7335_e7500_-e7710_e7730.pdf?rev=2e735fd18b544d46b36ee0e88353ef5c&sc_lang=en-us&hash=CC77B45817715F3ED3A8F3B1953450EB)
- Forward and reverse i5 adapter sequences.
- Forward and reverse i7 adapter sequences (modified based on i7 adapter sequence, nucleotide number, might vary between 6-8).

Note that the NNNNNN in the sequnces is the uniqe identifier for each sample/individual.
## Step 4: Sequence alignment, Sorting, and Cleaning
Run: `step4_align.sh` 

Software: `BWA` and `SAMtools`

Input: `fastq.gz` files, the output from trimmomatic, and reference genome (indexed with samtools faidx). Then the SAMtools uses the sam files generated from BWA to sort and clean the sequences.

Output: `.bam files`

Sorted files are smaller in size and faster to process; this is why the downstream tools require sorted files. The cleaning step removes aligned sequences with low scores.

Our lab already has the index for the reference genome in the same direcotry, so don't worry about making an index, but here is the script just in case:
```bash
#load BWA module 
module load BWA/0.7.17-GCCcore-12.2.0

# Path to reference genome
reference="/gpfs/gibbs/project/coughlan/hc858/test/Mguttatusvar_IM62_797_v3.0.fa"

# Create a directory for storing index files.
index_dir="ref_index_files"
mkdir -p $index_dir

# Change to the directory where the reference genome file is located.
cd $(dirname $reference)

# Run the BWA index command to generate the index and move the output to a new directory.
bwa index $(basename $reference)
mv *.bwt *.pac *.ann *.amb *.sa $index_dir
```
To check for the rate of mapping, you can run this for each bam
```bash
samtools flagstat 1A_1A_B6_sorted_clean.bam
```
Use the loop in the `calc_read_number_loop.sh` to see the read number in all generated .bam files
