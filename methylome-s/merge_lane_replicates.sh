#!/bin/bash -l
#SBATCH --job-name merge_lane_replicates
#SBATCH --requeue
#SBATCH --partition=general

########## Merge lane-replicate BAMs (e.g. L001 + L006) per sample #################
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

# run in analysis folder
cd analysis

########## Modules #################

module load samtools/1.13-gcc-11.3.0

# uncomment below (and comment out the module load above) to use a conda
# environment for samtools instead:
# conda activate $conda_enviro

########## Threads #################

THREADS=8

########## Set up sample #################

#get sample prefix from list
#use sed, -n suppression pattern space, then 'p' to print item number {SLURM_ARRAY_TASK_ID}
SAMPLE_PREFIX="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo "Processing sample: $SAMPLE_PREFIX"
echo "File suffix: $SUFFIX"
echo "Lane tags: $LANE1 + $LANE2"

# create output folder if it does not exist
mkdir -p $OUTPUT_DIR

# find_lane_bam <lane_tag>
# Matches by glob rather than exact concatenation, since the run ID and
# barcode between the sample prefix and lane tag differ per lane/run.
find_lane_bam() {
    local lane_tag=$1
    local matches=("${INPUT_DIR}/${SAMPLE_PREFIX}"_*_"${lane_tag}"*"${SUFFIX}")

    if [ ! -e "${matches[0]}" ]; then
        echo "ERROR: No BAM file found matching ${INPUT_DIR}/${SAMPLE_PREFIX}_*_${lane_tag}*${SUFFIX}"
        exit 1
    fi
    if [ "${#matches[@]}" -ne 1 ]; then
        echo "ERROR: Expected exactly one BAM file matching ${INPUT_DIR}/${SAMPLE_PREFIX}_*_${lane_tag}*${SUFFIX}, found ${#matches[@]}:"
        printf '%s\n' "${matches[@]}"
        exit 1
    fi
    echo "${matches[0]}"
}

LANE1_BAM="$(find_lane_bam "$LANE1")"
LANE2_BAM="$(find_lane_bam "$LANE2")"

echo "Lane 1 BAM: $LANE1_BAM"
echo "Lane 2 BAM: $LANE2_BAM"

OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_PREFIX}${SUFFIX}"

########## Run #################

echo "Starting samtools merge for $SAMPLE_PREFIX at $(date)"

# -f    overwrite output if it already exists
# -c    combine @RG headers with identical IDs (lane BAMs share the sample's RG)
# -@    number of additional threads to use
samtools merge -f -c -@ ${THREADS} "$OUTPUT_BAM" "$LANE1_BAM" "$LANE2_BAM"
samtools index -@ ${THREADS} "$OUTPUT_BAM"

# Check if output was created successfully
if [ -f "$OUTPUT_BAM" ] && [ -f "${OUTPUT_BAM}.bai" ]; then
    echo "SUCCESS: Output file created: $OUTPUT_BAM"
    read_count=$(samtools view -@ ${THREADS} -c "$OUTPUT_BAM")
    echo "Reads: $read_count"
else
    echo "ERROR: Output BAM or index was not created: $OUTPUT_BAM"
    exit 1
fi

echo "Completed samtools merge for $SAMPLE_PREFIX at $(date)"
