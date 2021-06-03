#!/bin/bash
#set -xe
set -xeuo pipefail

step=macs2

############ Change log ###########
# 26 May 2021 Misha Mangila
# Original code. Added help option
#
#
#
#

##################### Help ##########
Help()
{

echo "This script submits a batch job to PBS in order to run MACS2 for peak calling before or after deduplication."
echo
echo "USAGE:"
echo "bash macs_qsub.sh -i sample_list.txt -r reads_folder -A department -m memory
-c number_of_nodes -t walltime -s genome_size -p paired_end [-B(road)]"
echo
echo "For genome_size, input the effective genome size, either in scientific or dec
imal notation."
echo "Default value: 2.7e9."
echo
echo "For paired_end, input 'yes' for paired-end reads and 'no' for single-end reads."
echo "Default value: no"
echo
echo "Default memory usage is 20GB, default CPU usage is 6 CPUs and default walltime
 is 2 hours."
echo
echo "Add the tag -B to activate broad peak analysis and the tag -f if the samples have been filtered using the script 02b-filter-dups_qsub.sh."

}

#### Defaults ######
mem=20
nodes=1
walltime=2:00:00
genome_size=hs
paired_end=no
broad=no
filter=no

while getopts ":i:r:o:A:m:c:t:s:p:Bfh" option; do
   case $option in
      i) sample_list=${OPTARG};;
      r) reads_folder=${OPTARG};;
      o) output=${OPTARG};;
      A) account=${OPTARG};;
      m) mem=${OPTARG};;
      c) nodes=${OPTARG};;
      t) walltime=${OPTARG};;
      s) genome_size=${OPTARG};;
      p) paired_end=${OPTARG};;
      c) control=${OPTARG};;
      B) broad=yes;;
      f) filter=yes;;
      h) # display Help
         Help
         exit;;
     \?) # incorrect option
         echo "Error: Invalid option"
         exit;;
   esac
done


#define stepo in the pipeline - should be the same name as the script
step=macs2

############## Setup ###############

if [ "$1" == "-h" || "$1" == "--help" ]
then

echo $usage
echo
echo $size_explain
echo $PE_explain

exit 0
fi

errors=0

if [[ ! $sample_list ]]
then
errors=$(($errors + 1))
echo "ERROR: No sample list"
fi

if [[ ! $account ]]
then
errors=$(($errors + 10))
echo "ERROR: No account string"
fi

if [[ ! $reads_folder ]]
then
errors=$(($errors + 100))
echo "ERROR: No reads folder"
fi

if [[ $errors -gt 0 ]]
then
echo "Proper usage:"
echo
Help
exit $errors
fi

echo "Running MACS2"
cat $sample_list

number_of_samples=`wc -l $sample_list | awk '{print $1}'`

qsub_t=1
if [[ "$number_of_samples" -gt 1 ]]
then

qsub_t="1-${number_of_samples}"

fi

echo "argument to be passed to qsub -J is '$qsub_t'"

mkdir -p logs

timestamp=$(date +%Y%m%d-%H%M%S)
log_folder=logs/${timestamp}_${step}
mkdir $log_folder

scriptdir="$(dirname $(readlink -f $0))"
script_to_qsub=~${scriptdir}/${step}.sh
cat $script_to_qsub > ${log_folder}/script.log
cat $0 > ${log_folder}/qsub_runner.log

qsub -J $qsub_t \
-l walltime=${walltime},nodes=1:ppn=${nodes},mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-A ${account} \
-v LIST=${sample_list},paired_end=$paired_end,genome_size=$genome_size,broad=$broad,reads_folder=${reads_folder},filter=$filter
$script_to_qsub
