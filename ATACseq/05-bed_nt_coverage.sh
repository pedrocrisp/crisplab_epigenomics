#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N bigWigs_tiles
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
module load bedtools
# module load R/3.5.0-gnu
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p coverage_bed

########## Run #################

# bam_dir="trimmed_align_bowtie2"
# ID="UMR_McrBC_MF_90"
# chrom_sizes="/90days/uqpcrisp/tmp_refseqs/Nbenth/NbLab330.genome.chrom.sizes"
# reference_tile_file="/90days/uqpcrisp/tmp_refseqs/Nbenth/sites/NbLab330.genome_100bp_tiles.bed"
# CPM_filter_UMRs=0.5

# bedgraph nt resolution
# -pc Calculates coverage of intervals from left point of a pair reads to the right point.
bedtools genomecov -d -ibam $bam_dir/${ID}_sorted.bam -g $chrom_sizes |
awk -F$"\\t" 'BEGIN { OFS = FS } { print $1, $2-1, $2, $3 }' - |
sort -k1,1 -k2,2n - > coverage_bed/${ID}_nt_cov_sorted.bed

echo finished summarising
