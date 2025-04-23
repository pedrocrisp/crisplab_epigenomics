#!/bin/bash -l
#SBATCH --job-name featureCounts
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
# conda create --name featurecounts
# conda activate featurecounts
# conda install -c bioconda subread
conda activate $conda_enviro
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

outdir=analysis/${alignFolder}_featureCounts
mkdir -p ${outdir}

# check if single or paired end by looking for R2 file
# currently a read that overlaps multiple features at one location ie overlapping genes, is conuted for both genes (-O)
# reads mapping to multiple locations are excluded
if ([ ${format} == "SE" ])
then
echo single reads
featureCounts \
-F SAF \
-s $strand \
-O \
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
-O \
-a $reference \
-T 6 \
-o "$outdir/${ID}.counts" \
"${alignFolder}/${ID}.bam"
else
echo "ERROR: PE or SE not specified"
exit 1
fi

echo Done summarising
