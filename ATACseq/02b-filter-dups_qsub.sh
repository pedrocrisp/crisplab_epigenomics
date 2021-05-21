#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 02b-filter-dups_qsub.sh <sample_list.txt> <reads_folder> <paired_end> <walltime> <mem> <account_department>"

#define stepo in the pipeline - should be the same name as the script
step=02b-filter-dups

######### Setup ################
sample_list=$1
reads_folder=$2
paired_end=$3
walltime=$4
mem=$5
account_department=$6

if [ "$#" -lt "6" ]
then
echo $usage
exit -1
else
echo "initiating filtering"
cat $sample_list
echo genome reference is $bt2_genome
fi

#number of samples
# qsub -J only takes a range to hack for 1 samples create a second sample in the input list called "NULL" - should work?
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
qsub_t=1
else
qsub_t="1-${number_of_samples}"
fi
echo "argument to be passed to qsub -J is '$qsub_t'"

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
-l walltime=${walltime},nodes=1:ppn=2,mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-v LIST=${sample_list},reads_folder=$reads_folder,paired_end=$paired_end \
-A $account_department \
$script_to_qsub
