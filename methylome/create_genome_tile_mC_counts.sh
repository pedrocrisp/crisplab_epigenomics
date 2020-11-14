#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N create_genome_tile_mC_counts
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

module load R/3.5.0-gnu
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID
echo genome is $genome_prefix

mkdir -p sites
cd sites

########## Run MODULE 1 #################

# Run R module to creat 100bp tile bed file
R -f ~/gitrepos/crisplab_epigenomics/methylome/create_genome_tile_mC_counts_tiles.R \
--args $genome_prefix

##########              #################
########## Run MODULE 2 #################
##########              #################
# Run module to count number of CG/CHG/CHH sites in Each 100bp tile

module load bedtools/2.25.0

genome=$genome_prefix
echo $genome

bedtools getfasta -fi ../${genome}.fa \
-bed ${genome}_100pb_tiles_for_sites_calc.bed \
-fo ${genome}_100bp_tiles.fa

# these fasta records contain the 100 nt tile plus 2 nt either side for classifying sites bridging the tiles

####### counts sites
# make a fasta file of the motifs to search for
printf ">CG\nCG\n>CHG\nC[ATC]G\n>CHH\nC[ACT][ACT]" > contexts_motifs.fa

#cat contexts_motifs.fa
#>CG
#CG
#>CHG
#C[ATC]G
#>CHH
#C[ACT][ACT]

######### find motifs ##########

# seqkit locate finds the location of all CG/CHG and CHH sites
# awk removes any C sites identified in the first or last 2 bp flanks of each tile
# csvtk summarises per tile per context

time \
cat ${genome}_100bp_tiles.fa | \
seqkit locate -i -f contexts_motifs.fa | \
awk '(NR==1) ||
($1 ~ /:0-102/) && ($2 == "CG") && ($4 == "+") && ($5 >= 1) && ($6 <= 101) ||
($1 ~ /:0-102/) && ($2 == "CG") && ($4 == "-") && ($5 >= 1) && ($6 <= 100)  ||
($1 ~ /:0-102/) && ($2 == "CHG") && ($4 == "+") && ($5 >= 1) && ($6 <= 102) ||
($1 ~ /:0-102/) && ($2 == "CHG") && ($4 == "-") && ($5 >= 1) && ($6 <= 100)  ||
($1 ~ /:0-102/) && ($2 == "CHH") && ($4 == "+") && ($5 >= 1) && ($6 <= 104) ||
($1 ~ /:0-102/) && ($2 == "CHH") && ($4 == "-") && ($5 >= 1) && ($6 <= 100) ||
($1 !~ /:0-102/) && ($2 == "CG") && ($4 == "+") && ($5 >= 3) && ($6 <= 103) ||
($1 !~ /:0-102/) && ($2 == "CG") && ($4 == "-") && ($5 >= 2) && ($6 <= 102)  ||
($1 !~ /:0-102/) && ($2 == "CHG") && ($4 == "+") && ($5 >= 3) && ($6 <= 104) ||
($1 !~ /:0-102/) && ($2 == "CHG") && ($4 == "-") && ($5 >= 1) && ($6 <= 102)  ||
($1 !~ /:0-102/) && ($2 == "CHH") && ($4 == "+") && ($5 >= 3) && ($6 <= 104) ||
($1 !~ /:0-102/) && ($2 == "CHH") && ($4 == "-") && ($5 >= 1) && ($6 <= 102) ' | \
csvtk -t freq -k -f 1,2 > ${genome}_100bp_tiles_sites.tsv

####### counts sites inc N?
# make a fasta file of the motifs to search for
printf ">CG\nCG\n>CHG\nC[ATC]G\n>CHH\nC[ACT][ACT]\n>N\nN" > contexts_motifs_N.fa
#cat contexts_motifs_N.fa
#>CG
#CG
#>CHG
#C[ATC]G
#>CHH
#C[ACT][ACT]
#>N
#N

cat ${genome}_100bp_tiles.fa | \
seqkit locate -i -f contexts_motifs_N.fa | \
awk '(NR==1) ||
($1 ~ /:0-102/) && ($2 == "CG") && ($4 == "+") && ($5 >= 1) && ($6 <= 101) ||
($1 ~ /:0-102/) && ($2 == "CG") && ($4 == "-") && ($5 >= 1) && ($6 <= 100)  ||
($1 ~ /:0-102/) && ($2 == "CHG") && ($4 == "+") && ($5 >= 1) && ($6 <= 102) ||
($1 ~ /:0-102/) && ($2 == "CHG") && ($4 == "-") && ($5 >= 1) && ($6 <= 100)  ||
($1 ~ /:0-102/) && ($2 == "CHH") && ($4 == "+") && ($5 >= 1) && ($6 <= 104) ||
($1 ~ /:0-102/) && ($2 == "CHH") && ($4 == "-") && ($5 >= 1) && ($6 <= 100) ||
($1 !~ /:0-102/) && ($2 == "CG") && ($4 == "+") && ($5 >= 3) && ($6 <= 103) ||
($1 !~ /:0-102/) && ($2 == "CG") && ($4 == "-") && ($5 >= 2) && ($6 <= 102)  ||
($1 !~ /:0-102/) && ($2 == "CHG") && ($4 == "+") && ($5 >= 3) && ($6 <= 104) ||
($1 !~ /:0-102/) && ($2 == "CHG") && ($4 == "-") && ($5 >= 1) && ($6 <= 102)  ||
($1 !~ /:0-102/) && ($2 == "CHH") && ($4 == "+") && ($5 >= 3) && ($6 <= 104) ||
($1 !~ /:0-102/) && ($2 == "CHH") && ($4 == "-") && ($5 >= 1) && ($6 <= 102) ||
($1 ~ /:0-102/) && ($2 == "N") && ($6 <= 100) ||
($1 !~ /:0-102/) && ($2 == "N") && ($5 >= 3) && ($6 <= 102) ' | \
csvtk -t freq -k -f 1,2 > ${genome}_100bp_tiles_sites_Ns.tsv

##########              #################
########## Run MODULE 3 #################
##########              #################

# Run R module to count sites in each 100 bp tile
R -f ~/gitrepos/crisplab_epigenomics/methylome/create_genome_tile_mC_counts.R \
--args $genome_prefix

echo finished summarising
