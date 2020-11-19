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
mkdir -p bigWigs
mkdir -p tiles

########## Run #################

# bam_dir="trimmed_align_bowtie2"
# ID="UMR_McrBC_MF_90"
# chrc_sizes="/90days/uqpcrisp/tmp_refseqs/Nbenth/NbLab330.genome.chrom.sizes"
# tile_file="/90days/uqpcrisp/tmp_refseqs/Nbenth/sites/NbLab330.genome_100bp_tiles.bed"
# CPM_filter_UMRs=5

# bedgraph
bedtools genomecov -bg -ibam $bam_dir/${ID}_sorted.bam -g $chrc_sizes | sort -k1,1 -k2,2n - > bigWigs/${ID}_sorted.bedgraph
# bigwig
bedGraphToBigWig bigWigs/${ID}_sorted.bedgraph $chrc_sizes bigWigs/${ID}.bw

#### summarise coverage into 100bp tiles
# 5' coverage only
# the -a includes 0 counts, remove if not needed...
bedtools genomecov -bga -5 -ibam $bam_dir/${ID}_sorted.bam -g $chrc_sizes | sort -k1,1 -k2,2n - > bigWigs/${ID}_sorted_5prime.bedgraph

# use closest to overlap with 100bp tile then summarise with groupby
bedtools closest -t all -a bigWigs/${ID}_sorted_5prime.bedgraph -b $tile_file |
bedtools groupby -g 5,6,7 -c 4 > tiles/${ID}_100bp.bed

#### normalise to CPM
# norm factor (CPM) 4 decimal places
norm_factor=`awk 'BEGIN {sum1 += $4} END {print sum1/1000000}' tiles/${ID}_100bp.bed`
echo norm factor $norm_factor
#
awk -F$"\\t" -v norm_factor=$norm_factor 'BEGIN {OFS = FS} (NR>1){
  print $1, $2, $3, $4/norm_factor
}' "tiles/${ID}_100bp.bed" > tiles/${ID}_100bp_CPM.bed

#### filter
awk -v CPM_filter_UMRs=$CPM_filter_UMRs '$4 > CPM_filter_UMRs' tiles/${ID}_100bp_CPM.bed |
bedtools merge -c 4 -o sum -i - > tiles/${ID}_100bp_CPM_${CPM_filter_UMRs}_UMRs.bed

#### merge to call UMRs

echo finished summarising
