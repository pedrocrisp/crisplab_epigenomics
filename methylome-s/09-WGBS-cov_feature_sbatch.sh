#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 08-WGBS-cov_feature_qsub.sh <sample_list.txt> <walltime> <memory> <region_dir> <out folder> <CG_bigWig> <CHG_bigWig> <CHH_bigWig> <conda_enviro> <account_department>

This script is designed for determining the WGBS coverage over regions of interest
Where there is one reference WGBS bigwigs dataset but multiple region of interest files (eg genotypes or conditions or digests)
Each WGBS bigwig should be a single context.
It is run from the region project folder, eg if determining the WGBS coverage over a set of UMR-seq peaks, run in the UMR project folder.
"

#define stepo in the pipeline - should be the same name as the script
step=09-WGBS-cov_feature

######### Setup ################
sample_list=$1
walltime=$2
mem=$3
region_dir=$4
suffix=$5
out_dir=$6
CG_bigWig=$7
CHG_bigWig=$8
CHH_bigWig=$9
conda_enviro=${10}
account_department=${11}

if [ "$#" -lt "11" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for binning"
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
--cpus-per-task 2 \
--mem ${mem}gb \
-o ${log_folder}/${step}_o_%A_%a \
-e ${log_folder}/${step}_e_%A_%a \
--export LIST=${sample_list},region_dir=${region_dir},out_dir=${out_dir},suffix=${suffix},CG_bigWig=${CG_bigWig},CHG_bigWig=${CHG_bigWig},CHH_bigWig=${CHH_bigWig},conda_enviro=${conda_enviro} \
--account $account_department \
$script_to_sbatch

