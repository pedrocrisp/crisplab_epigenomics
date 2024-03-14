#!/bin/bash -l
#SBATCH --job-name fastqc
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

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo contrast is $ID
echo data folder is $path_to_data_files

########## Run #################

#Run R moudle to:
R -f ~/gitrepos/crisplab_epigenomics/methylome-s/11b-DSS_call_tile_DMRs.R \
--args ${ID} $DMR_contrasts_table_file $path_to_data_files

echo finished summarising
