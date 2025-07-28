# Segregation_Distortion_Mapping
follow: dx.doi.org/10.17504/protocols.io.bjnbkman for processing ddRAD data from raw fastqs to vcf (2021).
Make sure to request **i5** and **i7** fastq files from the sequencing center.
## Step 1: remove PCR duplicated from each library

This step removes PCR duplicate using the i5 molecular barcode.

python script: `step1_rmdup.py` (requires unzipped fastq files), runs in about 2 hours per library.
### Required inputs:
R1, R2, i5 and i7 fastq files
Rename files from sequencing center to match format of R1, R2, i5 and i7 fastq files:
files from YCGA will be supplied as:
I1 = I7 index
R1 = Sequencing Read 1
R2 = I5 index read
R3 = Sequencing Read 2

**To prefrom this step:** run `step1_rmdup.sh` script. This bash script will call the python script `step1_rmdup.py`.
This step will remove duplicates and output fastq.gz files for i5, i7, R1 and R2
The resulting files will have the same prefix, but will have the suffix .rmdup.1.fastq (forward reads) and .rmdup.2.fastq (reverse reads).

(skipping the step to flip the reads from the Fishman lab protocol, moving straight to next step to demultiplex samples in each library)

## Step 2: Sample Demultiplexing using Stacks
Requires a .txt file with sample barcodes of each well on the plate to demultiplex samples.
prepare a txt file listing sample IDs and corresponding barcodes. The barcode file is `A1_A1_barcodes`, note that the sample IDs will be different for rach library.

Important: The BestRAD protocol we used to construct library generates a unique “GG” at the beginning, so you have to add “GG” before you formal barcodes.

Important: for stacks, the input fastq files need to follow the naming format of fastq files as received from the sequencing center.


