#!/bin/bash
#SBATCH --job-name=1A_intrascaff
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=ycga
#SBATCH --time=45:00:00
#SBATCH --mail-type=ALL 
#SBATCH --mem=5G
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#note: This script cannot be threaded, so you only need 1 core.
#note: The run time of this script depends on the number of windows/markers and the number of individuals. For 463 individuals and ~3000 windows (i used 50 kb-sized windows), it took 40 to 45 hours. For the same number of windows but for 870 individuals, it took 90 hrs. Use the week partition instead of the ycga partition as the limit of the latter's limit is 48 hours.

module load SciPy-bundle/2020.11-foss-2020b-Python-2.7.18

python step11_hmm.intrascaff.R.py 1A_filtered_genotypes_err_rates.txt g.Genotypes.1A_2A_D8.bam.txt 1A_low_genotype_err_list.txt 1A_intrascaff_v1.txt
