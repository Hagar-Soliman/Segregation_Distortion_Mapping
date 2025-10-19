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

**⚠️Note:** This script DOES NOT filter the VCF and keeps multiallelic variants (i.e., does not filter to only keep biallelic sites). This is **important** as VCF filtering will take place downstream when calling ancestry. I also added two extra lines to create an index and a summary statistics file. The reason why this script does not do any filtering is that it can create a "raw" VCF that then can be filtered to different software depending on their requirements. 

## Step 6: VCF processing using 3 Python scripts
Run: `step6_VCF_prog1.py`, `step6_VCF_prog2.py`, and `step6_VCF_prog3.py` back to back.
- **step6_VCF_prog1.py:** This script processes a compressed VCF file containing SNP calls and filters variants based on quality and allele balance criteria. It reads a .vcf.gz file for a specified library (e.g., "1A"), extracts SNPs with a minimum mapping quality (MQ ≥ 20), and evaluates allele depth (AD) across samples. For each SNP, it calculates the number of samples with valid calls, average read depth, and average reference allele proportion. SNPs that pass thresholds for sample coverage (≥ 50 samples) and allele balance (between 0.2 and 0.8) are written to a slimmed-down VCF and two read depth summary files. It also tracks the number of SNPs per scaffold and the position of the last SNP, outputting this to a scaffold summary file. You can change the: Min_MQ_score, Min_lines_called, minQ (allele freq),maxQ(allele freq).
-  **step6_VCF_prog2.py:** This script processes SNPs scaffold-by-scaffold, evaluates read depth and allele balance (excluding two parental lines), and selects high-quality, well-distributed SNPs for further analysis—likely for RAD-tag based genotyping or linkage mapping. **⚠️Important Note:** Make sure to look at the .out file from the first Python program run and see which line contains the two parents and use this info to edit the second Python script (ln 62). Also, don't forget to change the number of individuals.
-  **step6_VCF_prog3.py:** This script processes a list of selected SNPs (SNPs.limited.txt) and extracts reference and alternate allele depths for each individual across those SNPs, organizing the output into scaffold-specific and sample-specific files. It will output a .txt file for all .bam files in their own directories, along with a .txt file for each chromosome (this will be output in the same directory where this script was ran from). Also, it should output an `all.samples.txt` file.

**⚠️Important Note:** that the third Python script outputs lots of data (.bam.txt file for each bam file in the original aligned directory) AND does not overwrite the generated files if it fails and is re-run. Instead, it will just keep adding lines to the pre-existing files. Thus, always make sure to remove any generated files if your script failed to run before re-running again.

## Step 7: Assign parental ancestry (R)
Run: `step7_assign_parental_ancestry.R` as a cluster job, as it will take a few hours to run. Use `step7_run_R_parental_ancestry.sh` to run. 

Before I start this step, I like to copy or move all the .bam.txt (including the parents and F1) files to a directory called, for example, 1A_refAlt and have the script and it's running bash file in the same directory. There is a loop to do so in the Loops file. Pia Schward edited this file so it can run parallel jobs to run faster.

This R script compares allele depths from two parental lines (IM62 and IMPO) to classify each SNP site as REF (IMPO), ALT(IM62), HET, or missing (NN), then uses this classification to polarize genotypes in a set of recombinant individuals, flipping allele counts where necessary to align with parental ancestry. For the polarization step, if parents are opposite (e.g., IMPO = REF, IM62 = ALT), it retains allele counts. - If parents are flipped (IMPO = ALT, IM62 = REF), swaps ref and alt counts to align with ancestry.

Output: `parental_ancestry.tt`, update the `.nam.txt' file to make sure they are poralized. and `all.scaffold_*txt.`

## step 8: Windows
Run: `step8_window.py`. Use `step8_window_run.sh`

This step averages across a window based on the number of SNPs/ number of reads/ some range. You can define heterozygote calls as het deviation (I changed it to 0.2 to make it stricter when calling heterozygous windows). This will output 3 files for each individual. 1) a Genotypes file, 2) a genostats file (total number of windows with AA, AB, NN genotypes), and 3) a windows file (which I think is a count for each site that was used to build the genotypes file).
**⚠️Important Note:** I hated how this script outputs files with no headers, so I added a few Python lines to add headers to the output files. However, downstream, these headers can cause an issue, which I tried to fix in the downstream script. For now, I added a # in front of those lines to just avoid the headache altogether. 

## Step 9: Second Ancestry loop to filter and output .g files (R)
This step will process the genotypes, genostats, and windows outputted from step 8. It will filter the genotypes files and also filter out windows with low depth. I manually filtered out bad individuals by concatenating all the genostats in one table (loop can be found in the Loops file) and excluded any individuals that:

**(1)** have 80% or of their windows missing or (NN)

**(2)** are 70% homozygous for the backcrossed parent or 70% heterozygous.

**(3)** are 70% homozygous for the other parent. Although BB windows are usually an AB window that was miscalled, but if it is at such high percentage, The indivdual is probably bad.  

I then made a new directory, copied all the files, and deleted all the individuals I excluded from the above criteria. I name this directory "1A_refAlt_filtered. After that, I ran the R script using this filtered directory. This is why the script did not pick up any bad individuals, just bad windows. 

Also, I ran this R script in the cluster's R line by line rather than a job like step 8. Update: When I made the window smaller (50kb) and had 963 indivduals, extracting the sites loops took almost 6 hrs. So, from now onward, I will submit the R script as a cluster job.

**⚠️Note:**  The Python script from step 8 will output files with a header if you decided to use the header line in the script, so make sure that header = TRUE if that was the case.

```r
#making sense of windows data:
#when het dev = 0.4 and windows= 100KB--> THIS IS WHAT I ENDED UP USING, BUT YOU MIGHT WANT TO PLAY WITH WINDOW SIZES

setwd("/home/hks25/palmer_scratch/VCF/1A_VCF/1A_refAlt_filtered")
path="/home/hks25/palmer_scratch/VCF/1A_VCF/1A_refAlt_filtered/"

files<-dir(path, pattern="Genostats.")
data<-data.frame(matrix(ncol=5,nrow=length(files)))
for(k in 1:length(files)){
  data[k,1:5]<-read.delim(files[k],header=TRUE)
}
names(data)<-c("indiv","AA","AB","BB","NN")
#check the data looks OK and they are the correct type
str(data)

#checking sites that suck:

genotype.files<-dir(path, pattern="Genotypes.")
countAA<-0
countAB<-0
countBB<-0
countNN<-0
test<-read.delim("Genotypes.1A_1A_A1.bam.txt", header=TRUE) #reads in random Genotypes file to get the site names
site.counts<-data.frame(matrix(nrow=2472,ncol=6)) #nrow=number of windows in your files aka the number of sites (get from the file in the line above)
window.names<-test[,2:3]
#This is the loop to start populating the site.counts table. The loops takes along time to run (~ 1hr for 416 indivduals)
View(site.counts) 
for(d in 1:2472){
  countAA<-0
  countAB<-0
  countBB<-0
  countNN<-0
  for(m in genotype.files){
    temp<-read.delim(m,header=TRUE)
    if(temp[d,4]=="AA"){
      countAA=countAA+1
    } else if(temp[d,4]=="AB"){
      countAB=countAB+1
    } else if(temp[d,4]=="BB"){
      countBB=countBB+1
    } else if(temp[d,4]=="NN"){
     countNN=countNN+1
    }
  }
  site.counts[d,3]<-countAA
  site.counts[d,4]<-countAB
  site.counts[d,5]<-countBB
  site.counts[d,6]<-countNN
}
site.counts[,1:2]<-window.names
names(site.counts)<-c("scaffold","pos","AA","AB","BB","NN")

#Write and save the table
write.table(site.counts, file="site.counts_hetdev=0.2+WS100K.txt",sep="\t",quote=FALSE,row.names=FALSE)
site.counts <- read.delim("/home/hks25/palmer_scratch/VCF/1A_VCF/1A_refAlt_filtered/site.counts_hetdev=0.2+WS100K.txt", header = TRUE)

#checking the number of sites with no coverage
count=0
for(i in 1:length(site.counts[,1])){
  if(site.counts[i,6]==540){ 
    count=count+1
  }
}

#found 363 windows/sites with no coverage in the all samples dataset

#number of individuals with no coverage
count=0
for(i in 1:length(data[,1])){
  if(data[i,5]==nrow(test)){
    count=count+1
  }
}

#there are 2 individuals with no coverage in the all samples dataset

#make some filtering decisions. Here, I've decided to keep any individual with <50 markers are removed (91 individuals) + cut sites if there's at least 345 individuals with missing data (~90%). These are probably too lax, but it's worth playing around with.

#For some reason, JKK's windows program duplicated a bunch of sites (aka certain scaffolds were repeated twice in the genotypes/windows file for each individual. Double check your output + remove the duplicate rows). In this set of files, duplicates for 100kb windows begin on 2170

#makes a list of sites to remove
site.list <- data.frame(matrix(nrow=0, ncol=2))
names(site.list) <- c("scaffold", "pos")

#Puprose:this code does is that it removes windows that are missing in 90% of the indivduals.
#Make sure to change count to 90% of total individuals, in this case, I have 416 invidiauls, thus, 0.9 X 416 = 375.
#change i to number of windows: here 2472, and change the threshold after if depedning on the sample size
for(i in 1:2472){
  if(site.counts[i, 6] > 375){
    site.list <- rbind(site.list, site.counts[i, 1:2])
  }
}
#note: HS found 543 windows in the list.
#write an output table of the sites list
write.table(site.list, file="1A_site_list0.2+WS100K_presentunder0.9.txt", sep="\t", quote=FALSE, row.names=FALSE)

#Purpose: Create a list indiv.list of individuals with low coverage to exclude from the analysis. Note: I have already done this (more or less) manually in excell by removing any inviduals if their NN windows is 80% of the total number of widnows.
#Condition: If the sum of AA, AB, and BB counts for an individual is less than 100, add that individual to indiv.list.

indiv.list <- c()
for(i in 1:length(data[, 1])){
  if((data[i, 2] + data[i, 3] + data[i, 4]) < 100){
    indiv.list <- append(indiv.list, data[i, 1])
  }
}

setwd("/home/hks25/palmer_scratch/VCF/1A_VCF/1A_refAlt_filtered/")
write.table(indiv.list, file="1A_indiv_list0.2+WS100K_presentunder100.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(indiv.list, file="1A_indiv_list.txt", sep="\t", quote=FALSE, row.names=FALSE)

## filter the Genotype.* files generated by JKK's windowing script:

#I'm resetting this from above because now I'm including the parents, but I'm calling the filelist something different just in case I need to go back.

setwd("/home/hks25/palmer_scratch/VCF/1A_VCF/1A_refAlt_filtered/")
path="/home/hks25/palmer_scratch/VCF/1A_VCF/1A_refAlt_filtered/"
output.path="/home/hks25/palmer_scratch/VCF/1A_VCF/1A_refAlt_filtered/1A_g_files/"
#made the indiv list based on the Genostats* files, so need to change those to read in the Genotypes* files. #

# Modify the filenames in indiv.list to include the prefix "Genotypes."
indiv.list <- paste0("Genotypes.", indiv.list, ".txt")

# Ensure that the filenames in indiv.list do not include paths, so they can match with geno.files
indiv.list <- basename(indiv.list)

#loops through Genotypes* files, only reads in individuals to keep (aka not on the list above), only records markers to keep (again, markers not on the list):

for(i in 1:length(genotype.files)){
  genotypes<-data.frame(matrix(nrow=0,ncol=4))
  if(!(genotype.files[i]%in%indiv.list)){
    temp<-read.delim(genotype.files[i],header=TRUE)
    for(j in 1:length(temp[1:2472,1])){
      f<-row.names(temp)
      if(!(f[j]%in%row.names(site.list))){
        genotypes<-rbind(genotypes,temp[j,])
      }
    }
    write.table(genotypes, file=paste(output.path,"g.",genotype.files[i],sep = ''), row.names=FALSE, col.names=FALSE, quote=FALSE, sep= "\t" )
  }
}
```
The above loop will produce the g. files needed for the next step. Notice that these files do not have a header. 

## Step 10 GOOGA: Galculate Genotype Error Rate
First, I moved all the .g files into a new directory called 1A_g_files. In this directory, have those scripts so they can run with no problems

Run: `step10_calc_genotype_err_rate.py` as a parallel job using this job list `step10_genotypes_list.txt`. This runs for each sample; use `step10_sample_list_loop.sh` to make the list.
```bash
module load dSQ
dsq --job-file step10_genotypes_list.txt --mem-per-cpu 4g -t 20:00 --mail-type ALL
#then run the sbatch code it will generate
```
Note: this Python script cannot run with barebone Python 2, this is why I am loading `SciPy-bundle/2020.11-foss-2020b-Python-2.7.18`

This Python script will generate a file for each .g file, and each file will contain one line. I used a loop called `combine_genotype_err_loop.sh` to combine them in a txt file called `1A_combined_genotypes_err.txt`. In this step, you will notice that there are errors in calculating the genotype error rate for some individuals. This happens, I believe, due to big gaps in some of the chromosomes that give that error. I have modified that Python script to prevent this error. But I still had some problematic individuals (just one or two) I exlcuded them and made a new list named `1A_filtered_genotypes_err_rates.txt`.

## Step 11 GOOGA: Calculate intra-scaffold recombinational fractions
Run: `step11_hmm.intrascaff.R.py` using `step11_intrascaff_rec_rates.sh`. This will take a few hours.

Input: `1A_filtered_genotypes_err_rates.txt`, a random `g.Genotypes.PlantID.bam.txt` to get the markers names, `1A_low_genotype_err_list.txt`, which is just a list of the individuals in the genotype error rates file (make sure it's the naked name with no `.bam.txt`. Finally, an empty `bad.marks.txt` 

The output file (1A_intrascaff_v1.txt) has your first linkage map! Column 1 is the chromosome, column 2 is the marker. In column 3, each value is the intrascaffold recombination rate between the marker for that row and the
marker below. The units on this are in probability of recombination (or cM/100). The maximum value here is 0.25 (this should not be a surprise because the max recombination rate is .5 (ie. unlinked) and there is a 50% chance that recombination occurs before the marker and a 50% chance it occurs after the marker (thus .25 each way). This may mean that a marker is either in the wrong location(in need of GOOGA!) or just poorly fits the genotype error algorithm and should have been discarded.

How do you determine which of these possibilities is correct? Well, we can discard these markers and see if there are still large recombination distances in the same spots. If it is the case, then our
chromosome needs reordering in Googa. If this fixes the problem, then the marker was simply bad.

Ok, let's get rid of the bad markers. So, highlight each of the crappy markers in yellow in Excel. There is likely either one maker with .25 or two makers. Remember, the recombination rate is calculated with the marker below it. So get rid of the marker below the .25 in case of a single marker with a .25. If there are two markers back to back with .25, then get rid of the second one. Compile the list of bad markers into your previously created tab-delimited text file named `bad.marks.txt`. You only need two columns, Chromosome and Marker, but do not put in a header line. Now we can rerun the program calculating intrascaffold error rates with the `bad.marks` version of the last program and see if we get rid of the 0.25 errors. It is normal for the last marker in the chromosome to have a large negative value, so don't delete that. 

## Step 12: GOOGA: Calculate Rntra-scaffold Recombinational Fractions After Removing Bad Markers
Run: `step12_hmm.intrascaff.R.badmarks.py` using `step12_intrascaff_rec_badmarks.sh`. This will take a few hours.

Input: `1A_filtered_genotypes_err_rates.txt`, a random genotype `g.Genotypes.1A_1A_D8.bam.txt`, a list of individuals "1A_low_genotype_err_list.txt`, and `bad.marks.txt`.

Output: `1A_intrascaff_v2.txt`

## Step 13: Estimate Genotypes at Each Locus for Each Individual
Before running this Python script, we will need to calculate the cumulative recombination fraction for each marker for each scaffold/chromosome. I used Excel to do so using the following instructions:
- First marker, cM = 0
- Second marker, cM = 100*D2
- Third maker, cM= 100*SUM($D$2:D3), then drag till the end of the chromosome
- Assuming your recombination fractions are in column D, with the first marker starting in cell D2. Start each chromosome over with this same code and then fill down to the end of the chromosome with the formula for the 3rd marker. Once done, copy and paste values only into a new sheet to avoid missing values. The final file should look like `1A_intrascaff_v3.txt`. The wrong format of this file can cause errors.

Run: `step13_run_genotype.pp.txt` that uses `step13_genotype.pp.py`. This should run as a parallel job on the cluster for each linkage group (i.e. Chr_01..14)

Output: This script outputs 14 files for each chromosome with posterior probabilities for each marker. `Chr_01..Chr_14.pp.txt`

I then use the following code to concatenate all of them. Make sure that all chromosomes are present after concatenating, because the length of this file (i.e. the number of rows) should equal the number of markers in `1A_intrascaff_v3.txt`
```bash
cat Chr_{01..14}.pp.txt > 1A_allmarker.txt
```

Then we move to R to make our genotype matrix and use rQTL to plot our data!!!!

The R script is called `1A_googa_processing.Rmd`. I usually do this on the cluster's R



