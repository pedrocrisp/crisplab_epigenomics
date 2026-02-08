#!/bin/bash
#SBATCH --job-name hisat2
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

module load hisat2/2.2.1--h87f3376_4/module

# *** updating samtools to gcc-11.3 (but v11 also being phased out in 2026 - test samtools/1.18-gcc-12.3.0) ***
module load samtools/1.13-gcc-11.3.0
# module load samtools/1.13-gcc-10.3.0
# module load samtools/1.18-gcc-12.3.0

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

outdir=analysis/${aligner}
mkdir -p ${outdir}

outdir_tmp=${TMPDIR}/${aligner}
mkdir -p ${outdir_tmp}

outsam="${outdir_tmp}/${ID}.sam"
outbam="${outdir}/${ID}.bam"
fastqcfolder=analysis/trimmed

# check how many satqs there are - assumes "fastq" suffix
fastqs="$(find $fastqcfolder -type f -name ${ID}*.fq.gz)"
# convert to array to count elements
fastqs_count=($fastqs)

# check if single or paired end by looking for R2 file
if ([ "${#fastqs_count[@]}" == 1 ] && [ ${aligner} == "hisat2" ] && [ ${strandedness} == "unstranded" ])
then
echo single reads unstranded
echo aligning with hisat2
hisat2 \
-x $index \
-p $threads \
-U $fastqs \
-S "$outsam"

elif ([ "${#fastqs_count[@]}" == 1 ] && [ ${aligner} == "hisat2" ])
then
echo single reads stranded
echo aligning with hisat2
hisat2 \
--rna-strandness $strandedness \
-x $index \
-p $threads \
-U $fastqs \
-S "$outsam"

elif ([ "${#fastqs_count[@]}" == 2 ] && [ ${aligner} == "hisat2" ] && [ ${strandedness} == "unstranded" ])
then
echo paired reads unstranded
echo aligning with hisat2
fq1="$(echo $fastqs |cut -d ' ' -f 1)"
fq2="$(echo $fastqs |cut -d ' ' -f 2)"
hisat2 \
-x $index \
-p $threads \
-1 $fq1 -2 $fq2 \
-S "$outsam"

elif ([ "${#fastqs_count[@]}" == 2 ] && [ ${aligner} == "hisat2" ])
then
echo paired reads stranded
echo aligning with hisat2
fq1="$(echo $fastqs |cut -d ' ' -f 1)"
fq2="$(echo $fastqs |cut -d ' ' -f 2)"
hisat2 \
--rna-strandness $strandedness \
-x $index \
-p $threads \
-1 $fq1 -2 $fq2 \
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
#Make an index of the sorted bam file (make a .csi index in case of big genomes - beware this could break other downstream tools...)
samtools index -c ${outbam}

# count total reads after filtering
echo "${ID} total alignments after MAPQ filter (should only include reads that map to one location)"
samtools view -c $outbam

rm $outsam
