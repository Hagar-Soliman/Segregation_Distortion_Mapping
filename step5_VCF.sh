#!/bin/bash
#SBATCH --job-name=1A_VCF
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL 
#SBATCH --mem=64G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load BCFtools/1.16-GCCcore-10.2.0
module load SAMtools/1.16

#first, compile a list of paths to all input bam files for variant calling into one bamfile list.txt

map=1A
REF_PATH=/gpfs/gibbs/project/coughlan/shared/genomes/ref/Mimulus_guttatus_var_IM62_v3.mainGenome.fasta
OUTPUT_DIR=/home/hks25/palmer_scratch/SNP_calling
bam_files=/home/hks25/palmer_scratch/SNP_calling/${map}_bam_list.txt


output_vcf=$OUTPUT_DIR/${map}.vcf.gz

echo "Calling SNPs"
bcftools mpileup -f "$REF_PATH" -a AD,FORMAT/DP -b "$bam_files" | \
bcftools call -m -v -Oz -o "$output_vcf"

# confirm that the VCF file has been successfully generated
if [ -f "$output_vcf" ]; then
    echo "VCF was generated：$output_vcf"
        # Index the VCF with tabix
        tabix -p vcf "$output_vcf"
        # Generate summary statistics for the VCF
        bcftools stats "$output_vcf" > "${OUTPUT_DIR}/${map}_vcf_stats.txt"
else
    echo "error：SNP calling failed，no VCF file was generated"
fi
