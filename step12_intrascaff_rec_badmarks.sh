#!/bin/bash
#SBATCH --job-name=1A_intrascaff_badmarks
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL 
#SBATCH --mem=16G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load SciPy-bundle/2020.11-foss-2020b-Python-2.7.18

python step12_hmm.intrascaff.R.badmarks.py 1A_filtered_genotypes_err_rates.txt g.Genotypes.1A_1A_D8.bam.txt 1A_low_genotype_err_list.txt 1A_intrascaff_v2.txt bad.marks.txt
