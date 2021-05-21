#!/bin/bash
#PBS -A UQ-SCI-SAFS
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

########## Modules #################

module load bowtie2/2.3.4.3
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

# For reads longer than about 50 bp Bowtie 2 is generally faster, more sensitive, and uses less memory than Bowtie 1.
# For relatively short reads (e.g. less than 50 bp) Bowtie 1 is sometimes faster and/or more sensitive.

# Bowtie 2 supports local alignment, which doesn't require reads to align end-to-end.
# Local alignments might be "trimmed" ("soft clipped") at one or both extremes in a way that optimizes alignment score.
# Bowtie 2 also supports end-to-end alignment which, like Bowtie 1, requires that the read align entirely.

#-x bowtie index
#--phred33
#--end-to-end dont trim reads to enable alignment
# --local allows soft clipping
#--mm memory-mapped I/O to load index allow multiple processes to use index
#-p number of threds to use
#-U fastq file path, and specifies reads are not paired
#-S file to write SAM alignemnts too (although default is stdout anyhow)
#-D and -R tell bowtie to try a little harder than normal to find alignments
#-L reduce substring length to 10 (default 22) as these are short reads
#-i reduce substring interval? more sensitive?
#-N max # mismatches in seed alignment; can be 0 or 1 (0)
#-D give up extending after <int> failed extends in a row (15)
# -k report N mapping locations
# --score-min L,0,0 sets formular for min aln score for alignment to be reported, eg readlength 36 * -0.6 + -0.6) = min 0; this means only report exact matches
# Bowtie 2 does not "find" alignments in any specific order,
# so for reads that have more than N distinct, valid alignments,
# Bowtie 2 does not guarantee that the N alignments reported are the best possible in terms of alignment score.
# Still, this mode can be effective and fast in situations where the user cares more about whether a read aligns
# (or aligns a certain number of times) than where exactly it originated.

if (( "${#fastqs_count[@]}" == 2 )); then

echo "paired reads"

# this is assume there is only one fq per sample - in R1 and R2 format for PE
# also assumes the suffix in "fq"
fq_1="${reads_folder}/${ID}_R1*.fq"
fq_2="${reads_folder}/${ID}_R2*.fq"

bowtie2 \
-x $bt2_genome \
--phred33 \
--local \
--very-sensitive-local \
-X 2000 \
-p $bt2_threads \
-1 $fq_1 \
-2 $fq_2 \
-S "$outsam"

else
echo "assuming single end"

bowtie2 \
-x $bt2_genome \
--phred33 \
--local \
--very-sensitive-local \
-X 2000 \
-p $bt2_threads \
-U $fastqs \
-S "$outsam"

fi

# count total alignments before after filtering
echo "${ID} total alignments before MAPQ filter (might include reads that map to multiple locations)"
samtools view -c ${outsam}

###### sort and index
# In Bowtie2, a read that aligns to more than 1 site equally well is never given higher than a MAPQ of 1
# (even if it aligns to only 2 sites equally well as discussed above).
# use mapq 10?
# if you want multialigning reads you will have lower this
samtools view -q 10 -b -@ $bt2_threads $outsam | samtools sort -m 8G -@ $bt2_threads -o $outbam

#Make an index of the sorted bam file
samtools index ${outbam}

# count total reads after filtering
echo "${ID} total alignments after MAPQ filter (should only include reads that map to one location)"
samtools view -c $outbam

#Delete the temporary sam.
rm -v ${outsam}
