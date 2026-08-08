#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash merge_lane_replicates_sbatch.sh <sample_list.txt> <file_suffix> <lane1_tag> <lane2_tag> <input_dir> <output_dir> <walltime> <memory> <conda_enviro> <account_department>

Example:
bash merge_lane_replicates_sbatch.sh samples.txt _sorted_MarkDup_pairs_clipOverlap.bam L001 L006 ./bams ./bams_merged 4:00:00 8 myenv my_department

Arguments:
  sample_list.txt     - File containing sample name prefixes (one per line), e.g. 'Creeper_D10'
                         (everything before the run ID in the raw filenames)
  file_suffix          - Common suffix shared by both lanes' BAM files (e.g. _sorted_MarkDup_pairs_clipOverlap.bam)
  lane1_tag            - First lane tag to merge (e.g. L001)
  lane2_tag            - Second lane tag to merge (e.g. L006)
  input_dir            - Directory containing the per-lane BAM files
  output_dir           - Directory for merged, indexed BAM files
  walltime             - Job walltime (HH:MM:SS format)
  memory               - Memory in GB per job
  conda_enviro         - Name of the conda environment to activate (must provide samtools)
  account_department   - SLURM account/department name

Merges each sample's two lane-replicate BAMs (e.g. L001 + L006, from the same
library sequenced across two runs/flow cells) into a single coordinate-sorted,
indexed BAM per sample using 'samtools merge'. One array task per sample.
Each lane's file is matched as '<sample_prefix>_*_<lane_tag><file_suffix>',
tolerating the run ID and barcode text between the sample prefix and lane tag.
Merged output is named '<sample_prefix><file_suffix>' so it matches the
'<sample_prefix>*<file_suffix>' glob used by downstream steps in this repo."

#define step in the pipeline - should be the same name as the script
step=merge_lane_replicates

######### Setup ################
sample_list=$1
file_suffix=$2
lane1_tag=$3
lane2_tag=$4
input_dir=$5
output_dir=$6
walltime=$7
mem=$8
conda_enviro=$9
account_department=${10}

if [ "$#" -lt "10" ]
then
echo "$usage"
exit -1
else
echo "Submitting samples listed in '$sample_list' for lane-replicate BAM merging"
echo "File suffix: '$file_suffix'"
echo "Lane tags: '$lane1_tag' + '$lane2_tag'"
echo "Input directory: '$input_dir'"
echo "Output directory: '$output_dir'"
echo ""
echo "Sample list:"
cat $sample_list
fi

#check that input directory exists
if [ ! -d "$input_dir" ]; then
    echo "ERROR: Input directory not found: $input_dir"
    exit -1
fi

#number of samples
# sbatch -J only takes a range; to handle 1 sample, just use '1'
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
sbatch_t=1
else
sbatch_t="1-${number_of_samples}"
fi
echo "Number of samples: $number_of_samples"
echo "Argument to be passed to sbatch --array is '$sbatch_t'"

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

#make output dir if it doesnt exist yet
mkdir -p $output_dir

#make logs folder, timestamped
log_folder=logs/${timestamp}_${step}
mkdir $log_folder

#script path and cat a record of what was run
script_to_sbatch=${scriptdir}/${step}.sh
cat $script_to_sbatch > ${log_folder}/script.log
cat $0 > ${log_folder}/sbatch_runner.log
cp $sample_list ${log_folder}/sample_list.txt

#submit sbatch and pass args
#-o and -e pass the file locations for std out/error
#--export additional variables to pass to the sbatch script including the array list and the dir structures
sbatch --array $sbatch_t \
-t ${walltime} \
-N 1 \
-n 1 \
--cpus-per-task 8 \
--mem ${mem}gb \
-o ${log_folder}/${step}_o_%A_%a \
-e ${log_folder}/${step}_e_%A_%a \
--export LIST=${sample_list},SUFFIX=${file_suffix},LANE1=${lane1_tag},LANE2=${lane2_tag},INPUT_DIR=${input_dir},OUTPUT_DIR=${output_dir},conda_enviro=${conda_enviro} \
--account $account_department \
$script_to_sbatch

echo ""
echo "Jobs submitted! Monitor with: squeue -u \$USER"
echo "Check logs in: $log_folder"
