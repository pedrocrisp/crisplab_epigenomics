#!/bin/bash -l
#SBATCH --job-name trim_galore_gz
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

module load fastqc/0.11.9-java-11
# cutadapt module causing issues - install v4.4 with conda and use that:
# conda create --name cutadapt python=3.7
# conda activate cutadapt
# conda install -c bioconda cutadapt
# conda install trim-galore
conda activate $conda_enviro

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {SLURM_ARRAY_TASK_ID} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

#make trimmed folder
trimmedfolder=analysis/trimmed
mkdir -p $trimmedfolder

fastqcfolder=analysis/fastqc
mkdir -p $fastqcfolder

# check if single or paired end or unzipped by looking for R2 file
if [ -e "reads/${ID}_R2_001.fastq.gz" ] || [ -e "reads/${ID}_R2.fastq.gz" ]; then

echo "paired reads"

#uncompress reads because trim_galore throws the error `gzip: stdout: Broken pipe` if I input .gz files
#gunzip reads/${ID}_R1_001.fastq.gz
#gunzip reads/${ID}_R2_001.fastq.gz

########## Run #################
# for swift libraries this trimms 20bp from the 5' end of the R2 read to remove the adaptase tail.
# Swift recommends symetrical trimming, so I trim from the R1 read too...

trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir $fastqcfolder" -o $trimmedfolder --paired reads/${ID}_R1*fastq.gz reads/${ID}_R2*fastq.gz

#compress original reads again
#gzip reads/${ID}_R1_001.fastq
#gzip reads/${ID}_R2_001.fastq

# check if single or paired end or unzipped by looking for R2 file
elif [ -e "reads/${ID}_R2_001.fastq" ] || [ -e "reads/${ID}_R2.fastq" ]; then

echo "paired reads uncompressed"

########## Run #################
# for swift libraries this trimms 20bp from the 5' end of the R2 read to remove the adaptase tail.
# Swift recommends symetrical trimming, so I trim from the R1 read too...

trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir $fastqcfolder" -o $trimmedfolder --paired reads/${ID}_R1*fastq reads/${ID}_R2*fastq

#compress original reads again
gzip reads/${ID}_R1_001.fastq
gzip reads/${ID}_R2_001.fastq

elif [ -e "reads/${ID}_R1_001.fastq.gz" ] || [ -e "reads/${ID}_R1.fastq.gz" ]; then
  # single end compressed
  ########## Run #################
  # for swift libraries this trimms 20bp from the 5' end of the R2 read to remove the adaptase tail.
  # Swift recommends symetrical trimming, so I trim from the R1 read too...

  trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir $fastqcfolder" -o $trimmedfolder reads/${ID}_R1*fastq.gz

else
echo "assuming single end uncompresed"

#uncompress reads because trim_galore throws the error `gzip: stdout: Broken pipe` if I input .gz files
#gunzip reads/${ID}_R1_001.fastq.gz

########## Run #################
# for swift libraries this trimms 20bp from the 5' end of the R2 read to remove the adaptase tail.
# Swift recommends symetrical trimming, so I trim from the R1 read too...

trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir $fastqcfolder" -o $trimmedfolder reads/${ID}_R1*fastq

#compress original reads again
gzip reads/${ID}_R1_001.fastq

fi

echo Done trimming
