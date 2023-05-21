#!/bin/bash
#set -xe
set -xeuo pipefail

#Peter Crisp
#2019-12-12
#Bash sbatch script for create_genome_tile_mC_counts.R

# this runs the 21- scripts but in a bacth submission loop per chr (hopefully)

usage="USAGE:
22-call-umrs-sbatch <sample_list> <reference_tile_file_folder> <chrom_sizes_folder> <chromomsome_list>
<coverage_filter_min>
<site_filter_min>
<MR_percent>
<UMR_percent>
<walltime> <memory> <cores> <account_department>
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
account_department=${12}

if [ "$#" -lt "12" ]
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

#sample_list_2=$(for SAMPLE in $(cat $sample_list); do
#echo ${SAMPLE}_${CHROMOSOME}
#done)

awk -F$"\\t" -v CHROMOSOME=$CHROMOSOME \
'BEGIN {OFS = FS} (NR=1){print $1"_"CHROMOSOME}' \
$sample_list > samples_${CHROMOSOME}.txt

#number of samples
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
sbatch_t=1
else
sbatch_t="1-${number_of_samples}"
fi
echo "argument to be passed to sbatch -t is '$sbatch_t'"

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
log_folder=logs/${timestamp}_${step}_${CHROMOSOME}
mkdir $log_folder

#script path and cat a record of what was run
script_to_sbatch=${scriptdir}/${step}.sh
cat $script_to_sbatch > ${log_folder}/script.log
cat $0 > ${log_folder}/sbatch_runner.log

#submit sbatch and pass args
#-o and -e pass the file locations for std out/error
#--export additional variables to pass to the sbatch script including the array list and the dir structures
sbatch --array $sbatch_t \
-t ${walltime} \
-N 1 \
-n 1 \
--cpus-per-task ${cores} \
--mem ${mem}gb \
-o ${log_folder}/${step}_o_%A_%a \
-e ${log_folder}/${step}_e_%A_%a \
--export LIST=samples_${CHROMOSOME}.txt,reference_100bp_tiles=$reference_100bp_tiles,chrom_sizes_path=$chrom_sizes_path,coverage_filter_min=$coverage_filter_min,site_filter_min=$site_filter_min,MR_percent=$MR_percent,UMR_percent=$UMR_percent \
--account $account_department \
$script_to_sbatch