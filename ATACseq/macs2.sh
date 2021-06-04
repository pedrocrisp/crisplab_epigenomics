#!/bin/bash
#PBS -N macs2
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

###### Modules ################
module load macs2

###### Set up directories #####
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

cd analysis
mkdir -p macs2

##### Errors ##################
if ([ $paired_end == "yes" ] || [ $paired_end == "no" ])
then
echo "Is paried end = " $paired_end
else
echo "Library type not specified correctly."
echo "Please indicate:"
echo "yes for PE"
echo "no for SE"
exit 1111
fi

if [ "$filter" == "yes" ]
then
input="-t ${reads_folder}_filtered/${ID}_sorted_NoDup_ProPairs.bam"
control_input="-c ${reads_folder}_filtered/${control}_sorted_NoDup_ProPairs.bam"
fi

if [ "$filter" == "no" ]
then
input="-t ${reads_folder}/${ID}_sorted.bam"
control_input="-c ${reads_folder}/${control}_sorted.bam"
fi

macs2_cmd="macs2 callpeak \
$input \
-g $genome_size \
-n ${ID}_test \
--outdir analysis/macs2"


if [ "$paired_end" == "yes" ]
then
macs2_cmd="$macs2_cmd -f BAMPE"
fi

if [ "$paired_end" == "no" ]
then
macs2_cmd="$macs2_cmd -f BAM"
fi

if [ "$broad" == "yes" ]
then
macs2_cmd="$macs2_cmd --broad"
fi

if [ -n "$control" ]
then
macs2_cmd="$macs2_cmd $control_input"
fi

$macs2_cmd

echo "MACS2 finished running"
exit 0
