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

**⚠️Note:** this Python script will not recognize your fasta files unless their name is ${library}_R1.fastaq. The prefix must be the library name, and the suffix must be R1, R2, i7, or i5. Otherwise, it won't input it.

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

Here are stacks flags:
```bash
  p: path to a directory of files.
  P,--paired: files contained within the directory are paired.
  I,--interleaved: specify that the paired-end reads are interleaved in single files.
  i: input file type, either 'fastq', 'gzfastq' (gzipped fastq), 'bam', or 'bustard' (default: guess, or gzfastq if unable to).
  b: path to a file containing barcodes for this run, omit to ignore any barcoding.
  o: path to output the processed files.
  f: path to the input file if processing single-end sequences.
  1: first input file in a set of paired-end sequences.
  2: second input file in a set of paired-end sequences.
  c,--clean: clean data, remove any read with an uncalled base ('N').
  q,--quality: discard reads with low quality (phred) scores.
  r,--rescue: rescue barcodes and RAD-Tags.
  t: truncate final read length to this value.
  D: capture discarded reads to a file.
  E: specify how quality scores are encoded, 'phred33' (Illumina 1.8+/Sanger, default) or 'phred64' (Illumina 1.3-1.5).
  w: set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).
  s: set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).
  y: output type, either 'fastq', 'gzfastq', 'fasta', or 'gzfasta' (default: match input type).

  Barcode options:
    --inline-null:   barcode is inline with sequence, occurs only on single-end read (default).
    --index-null:    barcode is provded in FASTQ header (Illumina i5 or i7 read).
    --null-index:    barcode is provded in FASTQ header (Illumina i7 read if both i5 and i7 read are provided).
    --inline-inline: barcode is inline with sequence, occurs on single and paired-end read.
    --index-index:   barcode is provded in FASTQ header (Illumina i5 and i7 reads).
    --inline-index:  barcode is inline with sequence on single-end read and occurs in FASTQ header (from either i5 or i7 read).
    --index-inline:  barcode occurs in FASTQ header (Illumina i5 or i7 read) and is inline with single-end sequence (for single-end data) on paired-end read (for paired-end data).

process_radtags 2.59
process_radtags -p in_dir [--paired [--interleaved]] [-b barcode_file] -o out_dir -e enz [-c] [-q] [-r] [-t len]
process_radtags -f in_file [-b barcode_file] -o out_dir -e enz [-c] [-q] [-r] [-t len]
process_radtags -1 pair_1 -2 pair_2 [-b barcode_file] -o out_dir -e enz [-c] [-q] [-r] [-t len]

  p: path to a directory of files.
  P,--paired: files contained within the directory are paired.
  I,--interleaved: specify that the paired-end reads are interleaved in single files.
  i: input file type, either 'fastq', 'gzfastq' (gzipped fastq), 'bam', or 'bustard' (default: guess, or gzfastq if unable to).
  b: path to a file containing barcodes for this run, omit to ignore any barcoding.
  o: path to output the processed files.
  f: path to the input file if processing single-end sequences.
  1: first input file in a set of paired-end sequences.
  2: second input file in a set of paired-end sequences.
  c,--clean: clean data, remove any read with an uncalled base ('N').
  q,--quality: discard reads with low quality (phred) scores.
  r,--rescue: rescue barcodes and RAD-Tags.
  t: truncate final read length to this value.
  D: capture discarded reads to a file.
  E: specify how quality scores are encoded, 'phred33' (Illumina 1.8+/Sanger, default) or 'phred64' (Illumina 1.3-1.5).
  w: set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).
  s: set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).
  y: output type, either 'fastq', 'gzfastq', 'fasta', or 'gzfasta' (default: match input type).

  Barcode options:
    --inline-null:   barcode is inline with sequence, occurs only on single-end read (default).
    --index-null:    barcode is provded in FASTQ header (Illumina i5 or i7 read).
    --null-index:    barcode is provded in FASTQ header (Illumina i7 read if both i5 and i7 read are provided).
    --inline-inline: barcode is inline with sequence, occurs on single and paired-end read.
    --index-index:   barcode is provded in FASTQ header (Illumina i5 and i7 reads).
    --inline-index:  barcode is inline with sequence on single-end read and occurs in FASTQ header (from either i5 or i7 read).
    --index-inline:  barcode occurs in FASTQ header (Illumina i5 or i7 read) and is inline with single-end sequence (for single-end data) on paired-end read (for paired-end data).

  Restriction enzyme options:
    -e <enz>, --renz-1 <enz>: provide the restriction enzyme used (cut site occurs on single-end read)
    --renz-2 <enz>: if a double digest was used, provide the second restriction enzyme used (cut site occurs on the paired-end read).
    Currently supported enzymes include:
      'aciI', 'ageI', 'aluI', 'apaLI', 'apeKI', 'apoI', 'aseI', 'bamHI', 
      'bbvCI', 'bfaI', 'bfuCI', 'bgIII', 'bsaHI', 'bspDI', 'bstYI', 'btgI', 
      'cac8I', 'claI', 'csp6I', 'ddeI', 'dpnII', 'eaeI', 'ecoRI', 'ecoRV', 
      'ecoT22I', 'haeIII', 'hinP1I', 'hindIII', 'hpaII', 'hpyCH4IV', 'kpnI', 'mluCI', 
      'mseI', 'mslI', 'mspI', 'ncoI', 'ndeI', 'ngoMIV', 'nheI', 'nlaIII', 
      'notI', 'nsiI', 'nspI', 'pacI', 'pspXI', 'pstI', 'rsaI', 'sacI', 
      'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'speI', 'sphI', 'taqI', 'xbaI', or 
      'xhoI'

  Protocol-specific options:
    --bestrad: library was generated using BestRAD, check for restriction enzyme on either read and potentially tranpose reads.

  Adapter options:
    --adapter-1 <sequence>: provide adaptor sequence that may occur on the single-end read for filtering.
    --adapter-2 <sequence>: provide adaptor sequence that may occur on the paired-read for filtering.
      --adapter-mm <mismatches>: number of mismatches allowed in the adapter sequence.

  Output options:
    --retain-header: retain unmodified FASTQ headers in the output.
    --merge: if no barcodes are specified, merge all input files into a single output file.

  Advanced options:
    --filter-illumina: discard reads that have been marked by Illumina's chastity/purity filter as failing.
    --disable-rad-check: disable checking if the RAD cut site is intact.
    --len-limit <limit>: specify a minimum sequence length (useful if your data has already been trimmed).
    --barcode-dist-1: the number of allowed mismatches when rescuing single-end barcodes (default 1).
    --barcode-dist-2: the number of allowed mismatches when rescuing paired-end barcodes (defaults to --barcode-dist-1).
```

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
## Step 4: Sequence alignment, Sorting, Cleaning, and Indexing
Run: `step4_align.sh` 

Software: `BWA` and `SAMtools`

Input: `fastq.gz` files, the output from trimmomatic, and reference genome (indexed with samtools faidx). Then the SAMtools uses the sam files generated from BWA to sort and clean the sequences.

Output: `.bam files`

Sorted files are smaller in size and faster to process; this is why the downstream tools require sorted files. The cleaning step removes aligned sequences with low scores.

Our lab already has the index for the reference genome in the same direcotry, so don't worry about making an index.

