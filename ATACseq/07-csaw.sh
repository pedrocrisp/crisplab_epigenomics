#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N csaw
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
module load bedtools
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis

########## Run #################

## Run R module to creat 100bp tile bed file
#R -f ~/gitrepos/crisplab_epigenomics/ATACseq/07-csaw.R \
#--args $ID $DMR_contrasts_table_file $path_to_data_files $blacklist $filter_FC $bin_size $window_spacing

# Annotation module
# annotate the DE-UMR list with proximity to genes and TEs

########## Annotation #################

# Strategy
# 1. Overlap with gene loci and get distance from loci. Genes version is using v5 of the gene space referene: includes protein coding (syntenic/non-syntenic), miRNA, tRNA and lncRNA

# bedtools closest
# -a FILE
# -b FILE1 FILE2
# -d report distance (direction irrelevant when comparing 2 DMR lists
# -t all
# -d distance reported as final column
# -D b  Report distance with respect to B. When B is on the - strand, “upstream” means A has a higher (start,stop). Want to know if DMR is 5'/3' of genes etc
# -mdb all report the closest feature from the b files NOT the closest feature in each file
# NOTE: -t Controlling how ties for “closest” are broken (default all reported - the list will grow)

# bedtools overlap with the gene annotation
projectFolder=${path_to_data_files}_CSAW
DE_UMR_bed_folder=${ID}_FCF${filter_FC}_bin${bin_size}_space${window_spacing}
DE_UMR_bed_file=${ID}_FCF${filter_FC}_bin${bin_size}_space${window_spacing}_differential_UMRs_CSAW_sig_metadata.bed

outFolder=${projectFolder}/${DE_UMR_bed_folder}/annotation
# gene annotation

bedtools closest \
-a ${projectFolder}/${DE_UMR_bed_folder}/${DE_UMR_bed_file} \
-b ${genome_annotation} \
-mdb all \
-t all \
-D b \
-g ${chromosome_sizes} \
> ${outFolder}/${ID}_Olap_gene_space_V.bed

# genome_annotation: ~/umn/refseqs/maize/genome_annotations/Zea_mays_AGPv4_36_fixed_introns_gene_ncRNA_synteny_miRbase_space_V_stranded_bed6.bed
# chromosome_sizes: ~/umn/refseqs/maize/genome_annotations/Zea_mays.AGPv4.dna.toplevel_sorted.chrom.sizes

# TE annotation
bedtools closest \
-a {projectFolder}/${DE_UMR_bed_folder}/${DE_UMR_bed_file}  \
-b ${TE_annotation} \
-mdb all \
-t all \
-D b \
-g ${chromosome_sizes} \
> ${outFolder}/${ID}_Olap_TE_order.bed

# TE_annotation: ~/umn/refseqs/maize/genome_annotations/Zea_mays.AGPv4.dna.toplevel_sorted.chrom.sizes

# summarise and plot annotation
R -f ~/gitrepos/crisplab_epigenomics/ATACseq/07-csaw-annotate-gene-proximity.R \
--args $ID $DMR_contrasts_table_file $path_to_data_files $blacklist $filter_FC $bin_size $window_spacing

echo finished summarising
