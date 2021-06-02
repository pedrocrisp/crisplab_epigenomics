#!/bin/bash
#PBS -N bowtie2
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

###### Modules ################
# modules
module load macs2

###### Set up directories #####
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

cd analysis
mkdir -p macs2

if [[ $paired_end != "yes" && $paired_end != "no" ]]
then

echo "Library type not specified correctly."
echo "Please indicate:"
echo "yes for PE"
echo "no for SE"
exit 1

fi

if [ "$paired_end" == "yes" ]
then

macs2 callpeak \
-t ${reads_folder}/${ID}_sorted_NoDup.bam \
-f BAMPE \
-g $genome_size \
-n ${ID}_test \
--broad \
--outdir $output

fi

if [ "$paired_end" == "no" ]
then
macs2 callpeak \
-t ${reads_folder}/${ID}_sorted_NoDup.bam \
-f BAM \
-g $genome_size \
-n ${ID}_macs2 \
--broad \
--outdir $output

fi

echo "MACS2 finished running"
exit 0
