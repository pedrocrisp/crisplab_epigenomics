#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 05-summarise_MethylDackel_qsub.sh <sample_list.txt> <genome.fa> <chromosome.sizes.file> <walltime> <memory> <ChrC_name> <cores>
for example:
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/05-summarise_methylation-WGBS_qsub.sh \
single_sample.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
12:00:00 \
60 \
no \
8
"

#define stepo in the pipeline - should be the same name as the script
step=05-summarise_MethylDackel

######### Setup ################
sample_list=$1
genome_reference=$2
chrom_sizes_file=$3
walltime=$4
mem=$5
ChrC_name=$6
cores=$7

if [ "$#" -lt "6" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for trimming"
cat $sample_list
echo genome reference is $genome_reference
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
-l walltime=${walltime},nodes=1:ppn=1,mem=${mem}gb \
-o ${log_folder}/${step}_o^array_index^ \
-e ${log_folder}/${step}_e^array_index^ \
-v LIST=${sample_list},genome_reference=$genome_reference,chrom_sizes_file=$chrom_sizes_file,make_subcontext=$make_subcontext,cores=$cores \
$script_to_qsub
