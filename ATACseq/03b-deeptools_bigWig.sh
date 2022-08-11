#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N deeptools_tiles
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
conda activate py3.7
module load Java/1.8.0_45
# picard needs R
module load R/3.5.0-gnu

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p ${bam_dir}_bigWigs_deeptools
mkdir -p ${bam_dir}_insert_metrics

########## Run #################
# make normalised bigWigs
# ising --extendReads as suggested for contiguous mapping data like ChIP to map whole of PE frag
bamCoverage \
--bam $bam_dir/${ID}_sorted.bam \
-o ${bam_dir}_bigWigs_deeptools/${ID}.bw \
--binSize 10 \
--normalizeUsing CPM \
--extendReads \
-p 2

# also make megabase coverage plots to look at any genomie-wide trimmed_align_bowtie2_bigWigs_deeptools
bamCoverage \
--bam $bam_dir/${ID}_sorted.bam \
-o ${bam_dir}_bigWigs_deeptools/${ID}.bedgraph \
--binSize 1000000 \
--normalizeUsing CPM \
--outFileFormat bedgraph \
-p 2

############ get frag sizes ##############
# get frag sizes to check level of tagmentation using picard
# path to call picard jar - eg  - "/home/uqpcrisp/software/picard.jar"

java -jar $path_to_picard CollectInsertSizeMetrics \
I=$bam_dir/${ID}_sorted.bam \
O=${bam_dir}_insert_metrics/${ID}_insert_size_metrics.txt \
H=${bam_dir}_insert_metrics/${ID}_insert_size_histogram.pdf \
M=0.5

echo finished summarising
