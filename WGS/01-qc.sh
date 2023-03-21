#!/bin/bash
#SBATCH --job-name fastqc
#SBATCH --requeue
#SBATCH --partition=general

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo SBATCH: working directory is $SLURM_SUBMIT_DIR
echo SBATCH: job identifier is $SLURM_JOBID
echo SBATCH: array_ID is ${SLURM_ARRAY_TASK_ID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to SLURM_SUBMIT_DIR
cd "$SLURM_SUBMIT_DIR"
echo working dir is now $PWD

########## Modules #################

module load fastqc/0.11.9-java-11

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

fastqcfolder=analysis/fastqc_raw
mkdir -p $fastqcfolder

# check how many satqs there are - assumes "fastq" suffix
fastqs="$(find ./reads -type f -name ${ID}*.fastq*)"
# convert to array to count elements
fastqs_count=($fastqs)

# check if single or paired end by looking for R2 file
if (( "${#fastqs_count[@]}" == 2 )); then

echo "paired reads"

########## Run #################
fastqc -o $fastqcfolder reads/${ID}_R1*.fastq.gz reads/${ID}_R2*.fastq.gz

else
echo "assuming single end"

########## Run #################

fastqc -o $fastqcfolder reads/${ID}_R1*.fastq.gz

fi

echo Done QC now you should run multiqc in output directory to summarise
