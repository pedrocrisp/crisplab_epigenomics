#!/bin/bash
#set -xe
set -xeuo pipefail

#Peter Crisp
#2019-12-12
#Bash qsub script for create_genome_tile_mC_counts.R

# this runs the 21- scripts but in a bacth submission loop per chr (hopefully)

usage="USAGE:
22-call-umrs-qsub <sample_list> <reference_tile_file_folder> <chrom_sizes_folder> <chromomsome_list>
<coverage_filter_min>
<site_filter_min>
<MR_percent>
<UMR_percent>
<walltime> <memory> <cores>
"

#define stepo in the pipeline - should be the same name as the script
step=21-call-umrs

######### Setup ################
sample_list=$1
reference_tile_file_folder=$2
chrom_sizes_folder=$3
chromomsome_list=$4
coverage_filter_min=$5
site_filter_min=$6
MR_percent=$7
UMR_percent=$8
walltime=$9
mem=${10}
cores=${11}

if [ "$#" -lt "11" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for analysis"
cat $sample_list
fi

# try to make a loop to submit a batch of batch job per chromosome...

for CHROMOSOME in $(cat $chromomsome_list); do

reference_100bp_tiles=${reference_tile_file_folder}/${CHROMOSOME}_100bp_tiles_zBased_sites_counts.txt
chrom_sizes_path=${chrom_sizes_folder}/${CHROMOSOME}.chrom.sizes

sample_list_2=$(for SAMPLE in $(cat $sample_list); do
echo ${SAMPLE}_${CHROMOSOME}
done)

#sample_list_2=$(awk -F$"_" -v CHROMOSOME=$CHROMOSOME \
#'BEGIN {OFS = FS} (NR=1){print $1, CHROMOSOME}' \
#$sample_list)

#number of samples
number_of_samples=`wc -l $sample_list_2 | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
qsub_t=1
else
qsub_t="1-${number_of_samples}"
fi
echo "argument to be passed to qsub -t is '$qsub_t'"

#find script to run, makes it file system agnostic
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi

########## Run #################

#make log and analysis folders
#make logs folder if it doesnt exist yet
mkdir -p logs

timestamp=$(date +%Y%m%d-%H%M%S)

#make analysis dir if it doesnt exist yet
analysis_dir=analysis
mkdir -p $analysis_dir

#make trimmgalore logs folder, timestamped
log_folder=logs/${timestamp}_${step}
mkdir $log_folder

#script path and cat a record of what was run
script_to_qsub=${scriptdir}/${step}.sh
cat $script_to_qsub > ${log_folder}/script.log
cat $0 > ${log_folder}/qsub_runner.log

#submit qsub and pass args
#-o and -e pass the file locations for std out/error
#-v additional variables to pass to the qsub script including the PBS_array list and the dir structures
qsub -J $qsub_t \
-l walltime=${walltime},nodes=1:ppn=${cores},mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-v LIST=${sample_list_2},reference_100bp_tiles=$reference_100bp_tiles,chrom_sizes_path=$chrom_sizes_path,coverage_filter_min=$coverage_filter_min,site_filter_min=$site_filter_min,MR_percent=$MR_percent,UMR_percent=$UMR_percent \
$script_to_qsub
done
