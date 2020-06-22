#!/bin/bash
#PBS -A UQ-SCI-SAFS
#PBS -N bowtie1
#PBS -r y
#PBS -m abe
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

module load bowtie/1.1.2
module load samtools/1.9

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID
#sample_dir="${reads_folder}/${ID}"

# enter analysis folder
cd analysis

# check how many satqs there are
fastqs="$(find $reads_folder -type f -name ${ID}*.fq*)"
# convert to array to count elements
fastqs_count=($fastqs)

#make adaligned folder bowtie2 (caution this will not fail if dir already exists)
outdir="${reads_folder}_align_bowtie2"
mkdir -p ${outdir}

# output structure
outsam="${outdir}/${ID}.sam"
tmpbam="${outdir}/${ID}.bam"
outbam="${outdir}/${ID}_sorted.bam"

########## Run #################

# -X 1000 \ max insert
# -m 1 \ only report unique mapping reads (max one valid alignment, exclude others)
# -v 2 \ allow up to 2 mismatches *CONSIDER INCREASEING IF READS ARE LONGER*
# --best \ bowtie guarantees the reported alignment(s) are the best in terms of the number of mismatches
# –strata \ Specifying --strata in addition to -a and --best causes bowtie to report only those alignments in the best alignment stratum
# -p $bt2_threads \ threads
# -1 $fq_1 \ forward
# -2 $fq_2 \ reverse
# -t \ report timings
# $bt1_genome \ index
# -S "$outsam" out sam

if (( "${#fastqs_count[@]}" == 2 )); then

echo "paired reads"

# this is assume there is only one fq per sample - in R1 and R2 format for PE
# also assumes the suffix in "fq"
fq_1="${reads_folder}/${ID}_R1*.fq"
fq_2="${reads_folder}/${ID}_R2*.fq"

bowtie \
-X 1000 \
-m 1 \
-v 2 \
--best \
–strata \
-p $bt2_threads \
-1 $fq_1 \
-2 $fq_2 \
-t \
$bt1_genome \
-S "$outsam"

else
echo "assuming single end"

bowtie \
-X 1000 \
-m 1 \
-v 2 \
--best \
–strata \
-p $bt2_threads \
-1 $fastqs \
-t \
$bt1_genome \
-S "$outsam"

fi


###### sort and index
# consider filtering: -q 10 is MAPQ >= 10; 1 in 10 chance mapping location is wrong for example
samtools view -b -@ $bt1_threads $outsam | samtools sort -m 8G -@ $bt1_threads -o $outbam

#Make an index of the sorted bam file
samtools index ${outbam}

#Delete the temporary sam.
rm -v ${outsam}
