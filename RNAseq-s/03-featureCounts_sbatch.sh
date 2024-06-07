#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 03-featureCounts_sbatch.sh <sample_list.txt> <walltime> <memory> <alignFolder> <strand> <format> <reference> <account_department>"

#define stepo in the pipeline - should be the same name as the script
step=03-featureCounts

######### Setup ################
sample_list=$1
walltime=$2
mem=$3
alignFolder=$4
strand=$5
format=$6
reference=$7
account_department=$8

if [ "$#" -lt "8" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for alignment"
cat $sample_list
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
echo "argument to be passed to sbatch -t is '$sbatch_t'"

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
--cpus-per-task 6 \
--mem ${mem}gb \
-o ${log_folder}/${step}_o_%A_%a \
-e ${log_folder}/${step}_e_%A_%a \
--export LIST=${sample_list},strand=$strand,format=$format,reference=$reference,alignFolder=$alignFolder \
--account $account_department \
$script_to_sbatch

