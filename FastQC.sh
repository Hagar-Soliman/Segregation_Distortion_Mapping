#!/bin/sh
#SBATCH --partition=scavenge
#SBATCH --time=1:00:00
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH --job-name=1A_1A_fastqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagar.soliman@yale.edu
#SBATCH --output="./fastqc_logs/slurm-%A_%a.out"
#SBATCH --error="./fastqc_logs/slurm-%A_%a.err"
#SBATCH --requeue
#SBATCH --array 0-47


#Specifies that the job will be requeued after a node failure.

#The default is that the job will not be requeued.
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory"=$SLURM_SUBMIT_DIR

library=1A_1B
OUTPUT_DIR="/home/hks25/palmer_scratch/libs/fastqc/"
INPUT_DIR="/home/hks25/palmer_scratch/libs/demultiplexed/${library}_demultiplexed"

mkdir -p "${OUTPUT_DIR}"
cd ${OUTPUT_DIR}

#generate a results directory

mkdir -p ./${library}_fastqc_results
mkdir -p ./${library}_fastqc_logs

#load required modules

module load FastQC/0.12.1-Java-11
ulimit -s unlimited
module load Java/11.0.

#search all the fastq files from the "data" directory and generate the array
files=($(find ${INPUT_DIR} -name "*.fq.gz" | sort))
file=${files[$SLURM_ARRAY_TASK_ID]}
fastqc -o ./${library}_fastqc_results/ ${file}
