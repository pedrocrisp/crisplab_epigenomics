#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 02-bowtie1_qsub.sh <sample_list.txt> <reads_folder> <bt1_threads> <bt1_genome.fa> <walltime> <mem>"

#define stepo in the pipeline - should be the same name as the script
step=02-bowtie1

######### Setup ################
sample_list=$1
reads_folder=$2
bt1_threads=$3
bt1_genome=$4
walltime=$5
mem=$6
if [ "$#" -lt "6" ]
then
echo $usage
exit -1
else
echo "initiating bowtie jobs on $reads_folder folder, bowtie can use $bt1_threads threads"
cat $sample_list
echo genome reference is $bt1_genome
fi

#number of samples
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
qsub_t="1-${number_of_samples}"
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
-l walltime=${walltime},nodes=1:ppn=8,mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-v LIST=${sample_list},reads_folder=$reads_folder,bt1_threads=$bt1_threads,bt1_genome=$bt1_genome \
$script_to_qsub
