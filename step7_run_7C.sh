#!/bin/bash
#SBATCH --job-name=7C_R_ancestry
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=16G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load R/4.2.0-foss-2020b

Rscript step7_7C_assign_parental_ancestry.R
