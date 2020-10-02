#!/bin/bash
#PBS -N 05-subread
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

module load samtools/1.10

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

outdir="${aligner}"
mkdir ${outdir}
outsam="${outdir}/${ID}.sam"
fastqcfolder=analysis/trimmed

# check how many satqs there are - assumes "fastq" suffix
fastqs="$(find ./fastqcfolder -type f -name ${ID}*.fastq.gz)"
# convert to array to count elements
fastqs_count=($fastqs)

# check if single or paired end by looking for R2 file
if ([ "${#fastqs_count[@]}" == 2 ] && [ ${aligner} == "subread-align" ])
then
echo paired reads
echo aligning with subread-align
$aligner -T 6 -t 0 -i $index --SAMoutput -r $fastqs -o "$outsam"
elif ([ "${#fastqs_count[@]}" == 2 ] && [ ${aligner} == "subread-align" ])
then
fq1="$(echo $fastqs |cut -d ' ' -f 1)"
fq2="$(echo $fastqs |cut -d ' ' -f 2)"
$aligner -T 6 -t 0 -i $index --SAMoutput -r ${fq1} -R ${fq2} -o "$outsam"
elif ([ "${#fastqs_count[@]}" == 1 ] && [ ${aligner} == "subjunc" ])
then
$aligner -T 6 -i $index --SAMoutput -r $fastqs -o "$outsam"
elif ([ "${#fastqs_count[@]}" == 2 ] && [ ${aligner} == "subjunc" ])
then
fq1="$(echo $fastqs |cut -d ' ' -f 1)"
fq2="$(echo $fastqs |cut -d ' ' -f 2)"
$aligner -T 6 -i $index --SAMoutput -r ${fq1} -R ${fq2} -o "$outsam"
else
echo "ERROR: not able to align multiple fq files per pair"
echo "fastqs:"
echo "${fastqs}"
exit 1
fi

echo Done aligning
