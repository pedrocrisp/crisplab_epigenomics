#!/bin/bash
#PBS -N 02-hisat2
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

outdir=analysis/${aligner}
mkdir -p ${outdir}
outsam="${outdir}/${ID}.sam"
outbam="${outdir}/${ID}.bam"
fastqcfolder=analysis/trimmed

# check how many satqs there are - assumes "fastq" suffix
fastqs="$(find $fastqcfolder -type f -name ${ID}*.fq.gz)"
# convert to array to count elements
fastqs_count=($fastqs)

# check if single or paired end by looking for R2 file
if ([ "${#fastqs_count[@]}" == 1 ] && [ ${aligner} == "hisat2" ])
then
echo single reads
echo aligning with hisat2
hisat2 \
--rna-strandness $strandedness \
-x $index \
-p $threads \
-U $fastqs \
-S "$outsam"

elif ([ "${#fastqs_count[@]}" == 2 ] && [ ${aligner} == "hisat2" ])
then
echo paired reads
echo aligning with hisat2
fq1="$(echo $fastqs |cut -d ' ' -f 1)"
fq2="$(echo $fastqs |cut -d ' ' -f 2)"
hisat2 \
--rna-strandness $strandedness \
-x $index \
-p $threads \
-1 fq1 -2 fq2 \
-S "$outsam"

else
echo "ERROR: not able to align multiple fq files per pair"
echo "fastqs:"
echo "${fastqs}"
exit 1
fi

echo Done aligning

echo "${ID} total alignments before MAPQ filter (might include reads that map to multiple locations)"
samtools view -c ${outsam}

# filter sam output and convert to bam
# sort
# -q 10 i think will remove multimapping, not sure anything >2 makes a difference until 60... - https://www.biostars.org/p/244106/
samtools view -bhu -q 10 -@ $threads $outsam | \
samtools sort -T ${ID} -m 2G -o $outbam -
#Make an index of the sorted bam file
samtools index ${outbam}

# count total reads after filtering
echo "${ID} total alignments after MAPQ filter (should only include reads that map to one location)"
samtools view -c $outbam

rm $outsam
