#!/bin/bash
#SBATCH --job-name fastqc
#SBATCH --requeue
#SBATCH --partition=general

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo "SBATCH: working directory is $(pwd)"
echo "SBATCH: job identifier is $SLURM_JOBID"
echo "SBATCH: array_ID is ${SLURM_ARRAY_TASK_ID}"
echo ------------------------------------------------------

########## Modules #################

module load fastqc/0.11.9-java-11

########## Set up directories #################

# Use 'sed' with -n option for suppressing pattern space and 'p' to print item number {PBS_ARRAYID} (e.g., 2 from {list})
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

fastqcfolder=analysis/fastqc_raw
mkdir -p $fastqcfolder

# Check how many fastq files there are - assumes "fastq" suffix
fastqs="$(find ${FASTQ_DIR} -type f -name ${ID}.fastq)"

# Convert to array to count elements
fastqs_count=($fastqs)

# Check if single or paired end by looking for the number of files that match the sample name
if (( "${#fastqs_count[@]}" == 2 )); then
  # Run FastQC on paired-end reads
  fastqc -o $fastqcfolder ${FASTQ_DIR}/${ID}_R1*.fastq ${FASTQ_DIR}/${ID}_R2*.fastq
  echo "Paired reads"  
else
  # Run FastQC on single-end reads
  fastqc -o $fastqcfolder ${FASTQ_DIR}/${ID}_R1*.fastq
  echo "Assuming single end"
fi

echo "Done QC. Now you should run multiqc in the output directory to summarize."
