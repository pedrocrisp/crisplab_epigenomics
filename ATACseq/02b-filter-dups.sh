#!/bin/bash
#PBS -A UQ-SCI-SAFS
#PBS -N filter
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

# requires java 1.8
module load Java/1.8.0_45
module load bamtools
module load samtools/1.9

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID
#sample_dir="${reads_folder}/${ID}"

# enter analysis folder
cd analysis
mkdir -p ${reads_folder}_filtered

if [ "$paired_end" == "yes" ]
then

########## Run #################

#mark duplicates
#requires sorted input - using samtools sort in bsmap step (co-ordinate sorted)
# if co-ordinate sorted then pairs where the mate is unmapped or has secondary alignment are not marked as duplicate
# ASSUME_SORTED=true because sorting performed with samtools but samtools doesnt seem to add this flag to the headder
java -jar /home/uqpcrisp/software/picard.jar NoDuplicates \
I=${reads_folder}/${ID}_sorted.bam \
O=${reads_folder}_filtered/${ID}_sorted_NoDup.bam \
METRICS_FILE=${reads_folder}_filtered/${ID}_NoDupMetrics.txt \
ASSUME_SORTED=true \
CREATE_INDEX=False \
REMOVE_DUPLICATES=true

# keep properly paired reads using bamtools package
# note that some reads marked as properly paired by bsmap actually map to different chromosomes
bamtools filter \
-isMapped true \
-isPaired true \
-isProperPair true \
-in ${reads_folder}_filtered/${ID}_sorted_NoDup.bam \
-out ${reads_folder}_filtered/${ID}_sorted_NoDup_pairs.bam

#index bam
# index
samtools index ${reads_folder}_filtered/${ID}_sorted_NoDup_ProPairs.bam

echo "done filtering"

elif [ "$paired_end" == "no" ]
then

# remove duplicate
java -jar /home/uqpcrisp/software/picard.jar NoDuplicates \
I=${reads_folder}/${ID}_sorted.bam \
O={reads_folder}_filtered/${ID}_sorted_NoDup.bam \
METRICS_FILE={reads_folder}_filtered/${ID}_NoDupMetrics.txt \
ASSUME_SORTED=true \
CREATE_INDEX=False \
REMOVE_DUPLICATES=true

# rename to keep script naming convention consistent
mv bsmapped_filtered/${ID}_sorted_NoDup.bam bsmapped_filtered/${ID}_sorted_NoDup_ProPairs.bam

# index
samtools index bsmapped_filtered/${ID}_sorted_NoDup_ProPairs.bam

echo "done filtering"

else

echo "library type not specified correctly, please indicate yes for PE or no for SE"

fi
