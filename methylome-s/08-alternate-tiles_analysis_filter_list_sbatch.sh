#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 08-tiles_analysis_filter_list_sbatch.sh <sample_list.txt> <data_folder> <tile_list> <filter_suffix> 
<minCHH sites> <minCHH coverage> <minCHG sites> <minCHG coverage> <minCG sites> <minCG coverage> <account_department>
for example:
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/10-tiles_analysis_filter_list_sbatch.sh \
samples.txt \
tiles \
analysis_02_tiles_SeqCap_meta_140_samples/chh_2x_cov_62_sample_tile_list.tsv (or FALSE) \
2_samp_MN \
4 2 \
4 2 \
4 2 \
account_department
"

#define stepo in the pipeline - should be the same name as the script
step=08-alternate-tiles_analysis_filter_list

######### Setup ################
walltime=$1
mem=$2
sample_list=$3
data_folder=$4
tile_list=$5
filter_suffix=$6
minCHHs=$7
minCHH_cov=$8
minCHGs=$9
minCHG_cov=${10}
minCGs=${11}
minCG_cov=${12}
account_department=${13}

if [ "$#" -lt "12" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for filtering"
cat $sample_list
fi

#number of samples
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
sbatch_t="1-2"
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
--cpus-per-task 4 \
--mem ${mem}gb \
-o ${log_folder}/${step}_o_%A_%a \
-e ${log_folder}/${step}_e_%A_%a \
--export LIST=${sample_list},data_folder=$data_folder,tile_list=$tile_list,filter_suffix=$filter_suffix,minCHHs=$minCHHs,minCHH_cov=$minCHH_cov,minCHGs=$minCHGs,minCHG_cov=$minCHG_cov,minCGs=$minCGs,minCG_cov=$minCG_cov \
--account $account_department \
$script_to_sbatch
