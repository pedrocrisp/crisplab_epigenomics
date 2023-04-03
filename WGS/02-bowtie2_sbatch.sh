#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 02-bowtie1_sbatch.sh <sample_list.txt> <reads_folder> <bt2_threads> <bt2_genome.fa> <q10filter> <walltime> <mem> <account_department>"

#define stepo in the pipeline - should be the same name as the script
step=02-bowtie2

######### Setup ################
sample_list=$1
reads_folder=$2
bt2_threads=$3
bt2_genome=$4
q10filter=$5
walltime=$6
mem=$7
account_department=$8

if [ "$#" -lt "8" ]
then
echo $usage
exit -1
else
echo "initiating bowtie jobs on $reads_folder folder, bowtie can use $bt2_threads threads"
cat $sample_list
echo genome reference is $bt2_genome
fi

#number of samples
# sbatch -J only takes a range to hack for 1 samples create a second sample in the input list called "NULL" - should work?
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
sbatch_t=1
else
sbatch_t="1-${number_of_samples}"
fi
echo "argument to be passed to sbatch -J is '$sbatch_t'"

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
--cpus-per-task ${bt2_threads} \
--mem ${mem}gb \
-o ${log_folder}/${step}_o_%A_%a \
-e ${log_folder}/${step}_e_%A_%a \
--export LIST=${sample_list},reads_folder=$reads_folder,bt2_threads=$bt2_threads,bt2_genome=$bt2_genome,q10filter=$q10filter \
--account $account_department \
$script_to_sbatch