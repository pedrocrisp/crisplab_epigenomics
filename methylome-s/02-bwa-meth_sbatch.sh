#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 02-bwa-meth_sbatch.sh <sample_list.txt> <genome.fa> <walltime> <memory> <threads> <conda_enviro> <account_department>"

#define step in the pipeline - should be the same name as the script
step=02-bwa-meth

######### Setup ################
sample_list=$1
genome_reference=$2
walltime=$3
mem=$4
cores=$5
conda_enviro=$6
account_department=$7

if [ "$#" -lt "7" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for mapping"
cat $sample_list
echo genome reference is $genome_reference
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
--cpus-per-task $cores \
--mem ${mem}gb \
-o ${log_folder}/${step}_o_%A_%a \
-e ${log_folder}/${step}_e_%A_%a \
--export LIST=${sample_list},genome_reference=$genome_reference,cores=$cores,conda_enviro=$conda_enviro \
--account $account_department \
$script_to_sbatch
