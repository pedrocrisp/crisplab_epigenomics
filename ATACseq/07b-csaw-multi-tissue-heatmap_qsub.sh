#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 07b-csaw-multi-tissue-heatmap_qsub.sh <regions_list.txt> <walltime> <memory> <sample_key_file> <path_to_data_files> <blacklist> <filter_FC> <bin_size> <window_spacing>
"

#define stepo in the pipeline - should be the same name as the script
step=07b-csaw-multi-tissue-heatmap

######### Setup ################
sample_list=$1
walltime=$2
mem=$3
sample_key_file=$4
path_to_data_files=$5
blacklist=$6
filter_FC=$7
bin_size=$8
window_spacing=$9

if [ "$#" -lt "9" ]
then
echo $usage
exit -1
else
echo "CSAW '$sample_list' "
cat $sample_list
fi

#number of samples
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
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
-l walltime=${walltime},nodes=1:ppn=1,mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-v LIST=${sample_list},sample_key_file=$sample_key_file,path_to_data_files=$path_to_data_files,blacklist=$blacklist,filter_FC=$filter_FC,bin_size=$bin_size,window_spacing=$window_spacing \
$script_to_qsub
