#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 02-hisat2_qsub.sh <account_department> <sample_list.txt> <walltime> <memory> <threads> <aligner> <genome_index> <strandedness>"

#define stepo in the pipeline - should be the same name as the script
step=02-hisat2

######### Setup ################
account_department=$1
sample_list=$2
walltime=$3
mem=$4
threads=$5
aligner=$6
index=$7
strandedness=$8

if [ "$#" -lt "8" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for alignment"
cat $sample_list
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
echo "argument to be passed to qsub -t is '$qsub_t'"

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
script_dir=~/gitrepos/crisplab_epigenomics/RNAseq
script_to_qsub=${script_dir}/${step}.sh
cat $script_to_qsub > ${log_folder}/script.log
cat $0 > ${log_folder}/qsub_runner.log

#submit qsub and pass args
#-o and -e pass the file locations for std out/error
#-v additional variables to pass to the qsub script including the PBS_array list and the dir structures
qsub -J $qsub_t \
-l walltime=${walltime},nodes=1:ppn=${threads},mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-v LIST=${sample_list},index=$index,aligner=$aligner,threads=$threads,strandedness=$strandedness \
-A $account_department \
$script_to_qsub
