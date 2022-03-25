#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N bwa-meth
#PBS -r y
#PBS -m abej
#PBS -M p.crisp@uq.edu.au

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: array_ID is ${PBS_ARRAY_INDEX}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

########## Modules #################
# activate conda with python 3 and the bwa-meth install
conda activate py3.7
module load bwa
module load samtools

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p bwa-meth

########## Run #################

# check if single or paired end by looking for R2 file
if [ -e "trimmed/${ID}_R2_001_val_2.fq.gz" ]; then

echo "paired reads"

bwameth.py \
--threads $cores \
--reference ${genome_reference} \
trimmed/${ID}_R1_001_val_1.fq.gz trimmed/${ID}_R2_001_val_2.fq.gz \
| samtools view -b - \
> bwa-meth/${ID}_tmp.bam

# sort
samtools sort bwa-meth/${ID}_tmp.bam > bwa-meth/${ID}.bam
#rm tmp unsorted file
rm bwa-meth/${ID}_tmp.bam
# get mapping stats
samtools flagstat bwa-meth/${ID}.bam > ${ID}_flagstat.txt
# print mapping stats (to log file) too
cat ${ID}_flagstat.txt

else
echo "assuming single end"

bwameth.py \
--threads $cores \
--reference ${genome_reference} \
trimmed/${ID}_R1_001_trimmed.fq.gz \
| samtools view -b - \
> bwa-meth/${ID}_tmp.bam

# sort
samtools sort bwa-meth/${ID}_tmp.bam > bwa-meth/${ID}.bam
#rm tmp unsorted file
rm bwa-meth/${ID}_tmp.bam
# get mapping stats
samtools flagstat bwa-meth/${ID}.bam > ${ID}_flagstat.txt
# print mapping stats (to log file) too
cat ${ID}_flagstat.txt

fi

echo Done mapping
