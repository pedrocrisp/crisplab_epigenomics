#!/bin/bash
#set -xe
set -xeuo pipefail

#Peter Crisp
#2019-12-12
#Bash sbatch script for create_genome_tile_mC_counts.R

usage="USAGE:
21-call-umrs-sbatch <sample_list> <reference_100bp_tiles> <chrom_sizes_path>
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
reference_100bp_tiles=$2
chrom_sizes_path=$3
coverage_filter_min=$4
site_filter_min=$5
MR_percent=$6
UMR_percent=$7
walltime=$8
mem=$9
cores=${10}
account_department=${11}

if [ "$#" -lt "11" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for analysis"
cat $sample_list
fi

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
log_folder=logs/${timestamp}_${step}
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