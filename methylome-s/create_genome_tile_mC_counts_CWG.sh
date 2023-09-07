#!/bin/bash -l
#SBATCH --job-name CWG
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

module load r/4.2.1-foss-2022a
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

genome=$ID
echo genome is $genome

mkdir -p sites
cd sites

########## Run MODULE 1 #################

# Run R module to creat 100bp tile bed file
R -f ~/gitrepos/crisplab_epigenomics/methylome-s/create_genome_tile_mC_counts_CWG_tiles.R \
--args $genome

##########              #################
########## Run MODULE 2 #################
##########              #################
# Run module to count number of CG/CTG/CCG sites in Each 100bp tile

module load bedtools

bedtools getfasta -fi ../${genome}.fa \
-bed ${genome}_100pb_tiles_for_sites_calc.bed \
-fo ${genome}_100bp_tiles.fa

# these fasta records contain the 100 nt tile plus 2 nt either side for classifying sites bridging the tiles

####### counts sites
# make a fasta file of the motifs to search for
printf ">CAG\nCAG\n>CTG\nCTG\n>CCG\nCCG" > contexts_motifs.fa

#cat contexts_motifs.fa
#>CAG
#CAG
#>CTG
#CTG#
#>CCG
#CCG

######### find motifs ##########

# seqkit locate finds the location of all CAG/CTG and CCG sites
# awk removes any C sites identified in the first or last 2 bp flanks of each tile
# csvtk summarises per tile per context

time \
cat ${genome}_100bp_tiles.fa | \
seqkit locate -i -r -f contexts_motifs.fa | \
awk '(NR==1) ||
($1 ~ /:0-102/) && ($2 == "CAG") && ($4 == "+") && ($5 >= 1) && ($6 <= 101) ||
($1 ~ /:0-102/) && ($2 == "CAG") && ($4 == "-") && ($5 >= 1) && ($6 <= 100)  ||
($1 ~ /:0-102/) && ($2 == "CTG") && ($4 == "+") && ($5 >= 1) && ($6 <= 102) ||
($1 ~ /:0-102/) && ($2 == "CTG") && ($4 == "-") && ($5 >= 1) && ($6 <= 100)  ||
($1 ~ /:0-102/) && ($2 == "CCG") && ($4 == "+") && ($5 >= 1) && ($6 <= 104) ||
($1 ~ /:0-102/) && ($2 == "CCG") && ($4 == "-") && ($5 >= 1) && ($6 <= 100) ||
($1 !~ /:0-102/) && ($2 == "CAG") && ($4 == "+") && ($5 >= 3) && ($6 <= 103) ||
($1 !~ /:0-102/) && ($2 == "CAG") && ($4 == "-") && ($5 >= 2) && ($6 <= 102)  ||
($1 !~ /:0-102/) && ($2 == "CTG") && ($4 == "+") && ($5 >= 3) && ($6 <= 104) ||
($1 !~ /:0-102/) && ($2 == "CTG") && ($4 == "-") && ($5 >= 1) && ($6 <= 102)  ||
($1 !~ /:0-102/) && ($2 == "CCG") && ($4 == "+") && ($5 >= 3) && ($6 <= 104) ||
($1 !~ /:0-102/) && ($2 == "CCG") && ($4 == "-") && ($5 >= 1) && ($6 <= 102) ' | \
csvtk -t freq -k -f 1,2 > ${genome}_100bp_tiles_sites.tsv


##########              #################
########## Run MODULE 3 #################
##########              #################

# Run R module to count sites in each 100 bp tile
R -f ~/gitrepos/crisplab_epigenomics/methylome-s/create_genome_tile_mC_counts_CWG.R \
--args $genome

echo finished summarising
