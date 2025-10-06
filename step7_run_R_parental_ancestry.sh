#!/bin/bash
#SBATCH --job-name=R_assign_ancestry
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=70G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


#runtime for 450 samples about 5min per sample, about 32 hours runtime , but note, that these files that duplicated entries present!!!

module load R/4.2.0-foss-2020b

Rscript assign_parental_ancestry.R
