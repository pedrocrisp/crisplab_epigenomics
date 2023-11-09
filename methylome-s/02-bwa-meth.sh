#!/bin/bash -l
#SBATCH --job-name bwa-meth
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
# activate conda with python 3 and the bwa-meth install
conda activate $conda_enviro
module load bwa/0.7.17-gcc-10.3.0
module load samtools/1.13-gcc-10.3.0

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {SLURM_ARRAY_TASK_ID} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p bwa-meth
mkdir -p bwa-meth-flagstat

########## Run #################

# check if single or paired end by looking for R2 file
if [ -e "trimmed/${ID}_R2_001_val_2.fq.gz" ]; then

echo "paired reads"

bwameth.py \
--threads $cores \
--reference ${genome_reference} \
trimmed/${ID}_R1_001_val_1.fq.gz trimmed/${ID}_R2_001_val_2.fq.gz \
| samtools view --threads $cores -b - \
> bwa-meth/${ID}_tmp.bam

# get mapping stats - use unfiltered bam
samtools flagstat -@ $cores bwa-meth/${ID}_tmp.bam > bwa-meth-flagstat/${ID}_flagstat.txt
# print mapping stats (to log file) too
cat bwa-meth-flagstat/${ID}_flagstat.txt

# make bam
# only keep mapped reads (-F 4)
# sort
samtools view --threads $cores -b -F 4 bwa-meth/${ID}_tmp.bam | \
samtools sort --threads $cores -m 3G - > bwa-meth/${ID}_sorted.bam

# check if single or paired end by looking for R2 file
elif [ -e "trimmed/${ID}_R2_val_2.fq.gz" ]; then

echo "paired reads"

bwameth.py \
--threads $cores \
--reference ${genome_reference} \
trimmed/${ID}_R1_val_1.fq.gz trimmed/${ID}_R2_val_2.fq.gz \
| samtools view --threads $cores -b - \
> bwa-meth/${ID}_tmp.bam

# get mapping stats - use unfiltered bam
samtools flagstat -@ $cores bwa-meth/${ID}_tmp.bam > bwa-meth-flagstat/${ID}_flagstat.txt
# print mapping stats (to log file) too
cat bwa-meth-flagstat/${ID}_flagstat.txt

# make bam
# only keep mapped reads (-F 4)
# sort
samtools view --threads $cores -b -F 4 bwa-meth/${ID}_tmp.bam | \
samtools sort --threads $cores -m 3G - > bwa-meth/${ID}_sorted.bam

else
echo "assuming single end"

bwameth.py \
--threads $cores \
--reference ${genome_reference} \
trimmed/${ID}_R1_001_trimmed.fq.gz \
| samtools view --threads $cores -b - \
> bwa-meth/${ID}_tmp.bam

# get mapping stats - use unfiltered bam
samtools flagstat -@ $cores bwa-meth/${ID}_tmp.bam > bwa-meth-flagstat/${ID}_flagstat.txt
# print mapping stats (to log file) too
cat bwa-meth-flagstat/${ID}_flagstat.txt

# make bam
# only keep mapped reads (-F 4)
# sort
samtools view --threads $cores -b -F 4 bwa-meth/${ID}_tmp.bam | \
samtools sort --threads $cores -m 3G - > bwa-meth/${ID}_sorted.bam

fi

#rm tmp unsorted sam file
# uncomment the below rm after script dev finished
rm bwa-meth/${ID}_tmp.bam

echo Done mapping
