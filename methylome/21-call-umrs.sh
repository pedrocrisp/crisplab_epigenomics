#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N 21-call-umrs
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
sample_to_crunch=$ID

echo sample being mapped is $ID

annotation_suffix=_mC_domains_II_cov_${coverage_filter_min}_sites_${site_filter_min}_MR_${MR_percent}_UMR_${UMR_percent}
annotation_suffix2=cov_${coverage_filter_min}_sites_${site_filter_min}_MR_${MR_percent}_UMR_${UMR_percent}

########## Run MODULE 1 DOMAINS #################

# Run R module to Call mC domains vII
R -f ~/gitrepos/crisplab_epigenomics/methylome/21-call-umrs-mod1.R \
--args $reference_100bp_tiles $ID $annotation_suffix $chrom_sizes_path $coverage_filter_min $site_filter_min $MR_percent $UMR_percent

# Notes
# Currently this script removes the organelles by using filter(!chr %in% c("Mt", "Pt")) - this may not catch all organelles
# post filtering may be required to remove organelles or other undesired contigs

########## Run MODULE 2 UMRs #################

# Merging UMTs and no data/no sites

# 1. Run bedtools closest on merged no data/no sites vs UMRs
# 2. Keep any no data tile that has distance = 1 and 2 closest UMRs
# 3. get those tiles and merge into UMRs
# 4. Calculate how many UMRs get consolidated and how the size categories totals change.
# 5. Filter 200 and under UMRs

########## Run Merge and sort UMRs #################
# Do I need to keep the domain column?

########## Merge #################
module load bedtools

#less ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs.bed
#echo ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs.bed
#
sort -k1,1 -k2,2n ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs.bed > ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_sorted.bed

bedtools merge -i ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_sorted.bed > ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_sorted_merge.bed

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_sorted.bed
# 1254448 * 100 / 1000000

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_sorted_merge.bed
# 321464

# less ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_sorted_merge.bed
# q

# Read into R and filter to add size column
R -f ~/gitrepos/crisplab_epigenomics/methylome/21-call-umrs-mod2.R \
--args $reference_100bp_tiles $ID $annotation_suffix $chrom_sizes_path $coverage_filter_min $site_filter_min $MR_percent $UMR_percent

########## Sort #################

sort -k1,1 -k2,2n ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_size.bed > \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_size_sorted.bed


########## Run Merge and sort NDs #################

sort -k1,1 -k2,2n ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs.bed > ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_sorted.bed

bedtools merge -i ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_sorted.bed > ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_sorted_merge.bed

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_sorted.bed
#3372705
#5030855-maize

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_sorted_merge.bed
#821221
#1520303-maize

### find NDs tile between UMTs

# I get a warning here due to the name chaecking function... I think its a false negative

bedtools closest \
-a ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_sorted_merge.bed \
-b ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_size_sorted.bed \
-t all \
-d \
-g $chrom_sizes_path \
> ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/NDs_Olap_UMTs.bed

# less ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/NDs_Olap_UMTs.bed

# Read into R and filter to NDs with distance 1 and two closest tiles. Adjacent tiles have a distance of 1
# Add a size colum
# and also a column for their co-ordinates
# remove black listed regions on second pass
R -f ~/gitrepos/crisplab_epigenomics/methylome/21-call-umrs-mod3.R \
--args $reference_100bp_tiles $ID $annotation_suffix $chrom_sizes_path $coverage_filter_min $site_filter_min $MR_percent $UMR_percent

### Merge NDs inbetween UMTs
# combine files and sort
cat \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_size_sorted.bed \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_between_UMTs.bed \
| sort -k1,1 -k2,2n > ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input.bed

# less ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input.bed
# q

# merge making a delimited list of cols 4 and 5 -
# col 4 this mark tiles with NDs in them and also retain metafeatures for subsequent filtering for distal tiles and
# col 5 this is the siz of the tile - I want to filter after merging to mark any tile that is more than 30% (?) NDs - probably have to work backwards to un merge these???

bedtools merge -i ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input.bed \
-c 4,5,6 -o collapse,collapse,collapse > \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_NDs.bed

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input.bed
# 329888 (maize 241801)

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_NDs.bed
# 231576 (maize 206976)

# (329888-231576)
# (329888-231576)/329888*100
# this reduced the numner of merged UMTs by 98312  or 30% (compared to 15.3% for maize)

# Read into R and filter UMTs that are too high % NDs
R -f ~/gitrepos/crisplab_epigenomics/methylome/21-call-umrs-mod4.R \
--args $reference_100bp_tiles $ID $annotation_suffix $chrom_sizes_path $coverage_filter_min $site_filter_min $MR_percent $UMR_percent

### Merge filtered NDs inbetween UMTs

# combine files and sort
cat \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_size_sorted.bed \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_between_UMTs_pct_filtered.bed \
| sort -k1,1 -k2,2n > ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input2.bed

# less ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input2.bed
# q
# merge making a delimited list of cols 4 and 5 -
# col 4 this mark tiles with NDs in them and also retain metafeatures for subsequent filtering for distal tiles and
# col 5 this is the siz of the tile - I want to filter after merging to mark any tile that is more than 30% (?) NDs - probably have to work backwards to un merge these???

bedtools merge -i ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input2.bed \
-c 4,5,6 -o collapse,collapse,collapse > \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_NDs_pct_filtered.bed

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/input2.bed
#329888 (maize 241801)

wc -l ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMTs_merge_NDs_pct_filtered.bed
#274411 (maize 206976)

#(329888-274411)
#(329888-274411)/329888*100
# this reduced the numner of merged UMTs by 55477 or 16.8% (compared to 15.3% for maize)

### Merge filtered NDs and tiles with data to get all regions included in analysis

# combine files and sort
cat \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_mC_domains_${annotation_suffix2}_tiles_with_data.bed \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_NDs_between_UMTs_pct_filtered_4col.bed \
| sort -k1,1 -k2,2n > ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_mC_domains_${annotation_suffix2}_tiles_with_data_inc_NDs_merged.bed

### UMR final list and Size annotation
R -f ~/gitrepos/crisplab_epigenomics/methylome/21-call-umrs-mod5.R \
--args $reference_100bp_tiles $ID $annotation_suffix $chrom_sizes_path $coverage_filter_min $site_filter_min $MR_percent $UMR_percent

############ put the key bed files in new output folder
# this step puts the final UMR list and the tiles with data bed files in a new folders
# the remaining files can probably just bed deleted ie delete the whole mC_UMT_annotation_beds folders

mkdir -p ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds_final

# UMRs
mv \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_UMRs_6col.bed \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds_final/${sample_to_crunch}_${annotation_suffix2}_UMRs_6col.bed

# tiles with data
mv \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_mC_domains_${annotation_suffix2}_tiles_with_data.bed \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds_final/

# tiles with data including the merged ND tiles
mv \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds/${sample_to_crunch}_mC_domains_${annotation_suffix2}_tiles_with_data_inc_NDs_merged.bed \
${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds_final/

# delete those remaining files
# if you ever want to look at the intermediate files again comment out the next line
rm -rv ${sample_to_crunch}${annotation_suffix}/mC_UMT_annotation_beds

echo finished summarising
