#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 11b-DSS_call_tile_DMRs_sbatch.sh <contrasts_list.txt> <DMR_contrasts_table_file> <path_to_data_files> <account_department>
for example:
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/11b-DSS_call_tile_DMRs_sbatch.sh \
DMR_tests_combos_all_list.txt \
DMR_tests_combos_all_table.tsv \
analysis/tiles_filtered_4C_2x \
a_crisp
"

#define stepo in the pipeline - should be the same name as the script
step=11b-DSS_call_tile_DMRs

######### Setup ################
walltime=$1
mem=$2
sample_list=$3
DMR_contrasts_table_file=$4
path_to_data_files=$5
account_department=$6

if [ "$#" -lt "6" ]
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
--export LIST=${sample_list},DMR_contrasts_table_file=$DMR_contrasts_table_file,path_to_data_files=$path_to_data_files \
--account $account_department \
$script_to_sbatch