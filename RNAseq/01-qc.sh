#!/bin/bash
#PBS -N fastqc
#PBS -r y
#PBS -m abej

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

module load fastqc/0.11.4

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

fastqcfolder=analysis/fastqc
mkdir -p $fastqcfolder

# check how many satqs there are - assumes "fastq" suffix
fastqs="$(find ./reads -type f -name ${ID}*.fastq*)"
# convert to array to count elements
fastqs_count=($fastqs)

# check if single or paired end by looking for R2 file
if (( "${#fastqs_count[@]}" == 2 )); then

echo "paired reads"

########## Run #################
fastqc -o $fastqcfolder $fastqs reads/${ID}*1.fastq.gz reads/${ID}*2.fastq.gz

else
echo "assuming single end"

########## Run #################

fastqc -o r$fastqcfolder $fastqs reads/${ID}*1.fastq.gz

fi

echo Done QC now you should run multiqc in output directory to summarise
