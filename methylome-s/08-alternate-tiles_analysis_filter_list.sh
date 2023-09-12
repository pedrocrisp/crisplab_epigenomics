#!/bin/bash -l
#SBATCH --job-name tile_filter
#SBATCH --requeue
#SBATCH --partition=general

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo SBATCH: working directory is $SLURM_SUBMIT_DIR
echo SBATCH: job identifier is $SLURM_JOBID
echo SBATCH: array_ID is ${SLURM_ARRAY_TASK_ID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to SLURM_SUBMIT_DIR
cd "$SLURM_SUBMIT_DIR"
echo working dir is now $PWD

########## Modules #################
module load r/4.2.1-foss-2022a

########## Set up dirs #################
output_folder=${data_folder}_filtered_${filter_suffix}
mkdir -p $output_folder

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {SLURM_ARRAY_TASK_ID} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID
echo data folder is $data_folder
echo tile list is $tile_list

########## Run #################
#Run R moudle to:
# 1. make CHH coverage files for summary stats analysis
R -f ~/gitrepos/crisplab_epigenomics/methylome/08-alternate-tiles_analysis_filter_list.R \
--args ${ID} $data_folder $tile_list $output_folder $minCHHs $minCHH_cov $minCHGs $minCHG_cov $minCGs $minCG_cov

echo finished summarising
