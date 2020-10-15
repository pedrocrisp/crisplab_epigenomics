#!/bin/bash
#PBS -N 03-featureCounts
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


########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

outdir=analysis/featureCounts
mkdir -p ${outdir}

# check if single or paired end by looking for R2 file
if ([ ${format} == "SE" ])
then
echo single reads
featureCounts \
-F SAF \
-s $strand \
-a $reference \
-T 6 \
-o "$outdir/${ID}.counts" \
"${alignFolder}/${ID}.bam"
elif ([ ${format} == "PE" ])
then
echo paired reads
featureCounts \
-F SAF \
-p \
-C \
-s $strand \
-a $reference \
-T 6 \
-o "$outdir/${ID}.counts" \
"${alignFolder}/${ID}.bam"
else
echo "ERROR: PE or SE not specified"
exit 1
fi

echo Done summarising
