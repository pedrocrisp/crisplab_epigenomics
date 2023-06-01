#!/bin/bash -l
#SBATCH --job-name bsmap
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

# bsmap requires samtools < 1.0.0 (note: now output is sam and running samtools manually on my own)
# module load samtools/0.1.18 # no longer an available module
PATH=~/software/bsmap-2.74/samtools:$PATH

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p bsmapped

########## Run #################

# align adapter trimmed datasets to B73 genome
        # -r 0: [0,1]   how to report repeat hits, 0=none(unique hit/pair only); 1=random one, default:1
        # -v 5: allow 5 mismatches (could also use -v 0.05 = 5% of read length)
        # -p 8: 8 threads/cores
        # -q 20: trim to q20

# check if single or paired end by looking for R2 file
if [ -e "trimmed/${ID}_R2_001_val_2.fq.gz" ]; then

echo "paired reads"

bsmap \
-a trimmed/${ID}_R1_001_val_1.fq.gz \
-b trimmed/${ID}_R2_001_val_2.fq.gz \
-d ${genome_reference} \
-o bsmapped/${ID}.sam \
-v 5 \
-r 0 \
-p $cores \
-q 20 \
-A $adapter_seq

elif [ -e "trimmed/${ID}_R2_val_2.fq.gz" ]; then

echo "paired reads"

bsmap \
-a trimmed/${ID}_R1_val_1.fq.gz \
-b trimmed/${ID}_R2_val_2.fq.gz \
-d ${genome_reference} \
-o bsmapped/${ID}.sam \
-v 5 \
-r 0 \
-p $cores \
-q 20 \
-A $adapter_seq


else
echo "assuming single end"

bsmap \
-a trimmed/${ID}_R1_001_trimmed.fq.gz \
-d ${genome_reference} \
-o bsmapped/${ID}.sam \
-v 5 \
-r 0 \
-p $cores \
-q 20 \
-A $adapter_seq

fi

echo Done mapping
