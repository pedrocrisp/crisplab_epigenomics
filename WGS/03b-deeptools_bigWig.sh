#!/bin/bash -l
#SBATCH --job-name bowtie2
#SBATCH --requeue
#SBATCH --partition=general

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo SBATCH: working directory is $SLURM_SUBMIT_DIR
echo SBATCH: job identifier is $SLURM_JOBID
echo SBATCH: array_ID is ${SLURM_ARRAY_TASK_ID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to SLURM_SUBMIT_DIR
cd "$SLURM_SUBMIT_DIR"
echo working dir is now $PWD

########## Modules #################
conda activate py3.7
module load java/11.0.16
# picard needs R
module load r/4.2.1-foss-2021a

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

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

java -jar ~/software/picard.jar CollectInsertSizeMetrics \
I=$bam_dir/${ID}_sorted.bam \
O=${bam_dir}_insert_metrics/${ID}_insert_size_metrics.txt \
H=${bam_dir}_insert_metrics/${ID}_insert_size_histogram.pdf \
M=0.5

echo finished summarising
