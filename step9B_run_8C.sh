#!/bin/bash
#SBATCH --job-name=8C_step9B_50kb
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=6G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load R/4.2.0-foss-2020b

Rscript step9B_8C_googa_prep.R
