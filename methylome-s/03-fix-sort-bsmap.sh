#!/bin/bash -l
#SBATCH --job-name fix_sort
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

#bsmap requires samtools < 1.0.0
# module load samtools/0.1.18 # no longer an available module
PATH=~/software/bsmap-2.74/samtools:$PATH

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

cd analysis/bsmapped

########## Run #################

# make bam
samtools view -bS ${ID}.sam > ${ID}.bam

# sort by read name (needed for fixsam)
samtools sort -n ${ID}.bam ${ID}_nameSrt

# fix mate pairs
samtools fixmate ${ID}_nameSrt.bam ${ID}_nameSrt_fixed.bam

# co-ordinate sort
samtools sort ${ID}_nameSrt_fixed.bam ${ID}_sorted

# index
samtools index ${ID}_sorted.bam

# remove intermediate files
rm ${ID}.sam ${ID}_nameSrt.bam ${ID}_nameSrt_fixed.bam ${ID}.bam

echo Done fixing
