#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 03-bed_bigWig_tiles_qsub.sh <sample_list.txt> <walltime> <memory> <bam_dir> <account_department>
"

#define stepo in the pipeline - should be the same name as the script
step=03b-deeptools_bigWig

######### Setup ################
sample_list=$1
walltime=$2
mem=$3
bam_dir=$4
account_department=$5

if [ "$#" -lt "5" ]
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
-l walltime=${walltime},nodes=1:ppn=2,mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-v LIST=${sample_list},bam_dir=${bam_dir} \
-A $account_department \
$script_to_qsub
