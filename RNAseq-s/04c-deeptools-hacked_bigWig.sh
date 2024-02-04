#!/bin/bash -l
#SBATCH --job-name bigWig
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
# dependencies
# conda env with deeptools-hacked and mosdepth installed
# conda create deeptools-hacked python=3.7
# conda activate deeptools-hacked
# conda install -c bioconda deeptools=3.5.1
# in "writeBedGraph.py"
# eg found in: ~/miniconda3/pkgs/deeptools-3.5.1-pyhdfd78af_1/site-packages/deeptools
# uncomment:
#  # uncomment these lines if fixed step bedgraph is required
#            if not np.isnan(value):
#                writeStart = start + tileIndex * self.binLength
#                writeEnd  =  min(writeStart + self.binLength, end)
#                _file.write(line_string.format(chrom, writeStart,
#                                               writeEnd, value))
#            continue
#

conda activate $conda_enviro

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p ${bam_dir}_bigWigs_deeptools
mkdir -p ${bam_dir}_insert_metrics
mkdir -p ${bam_dir}_100bp_tiles
mkdir -p ${bam_dir}_big_tiles

# NOTE: do we need to clip overlapping reads?
# compare to bams and see what the script actually does

########## Run #################
# make normalised bigWigs
# using --extendReads as suggested for contiguous mapping data like ChIP to map whole of PE frag
# try making 1bp bins for higher res coverage calculation 
# if we are going to use the .bw for pulling coverage, this may increase precision eg in correlations between samples?

## should we skip zeros when using bamCoverage???
## skipping zeros would make files much smaller
## skipping zeros might impact coverage comparsions between samples depending on tool used for making comparisons...

if [ "$highres" == "highres" ]
then
echo making high-res 1bp resolution bigwig files **caution** these are big files
# 1 bp bins for higher res
# use for measuring and comparing coverage
bamCoverage \
--bam $bam_dir/${ID}.bam \
-o ${bam_dir}_bigWigs_deeptools/${ID}.1.bw \
--binSize 1 \
--normalizeUsing CPM \
--extendReads \
-p 2
fi

# 10 bp bins for smaller files
# 10x smaller than the 1bp files
# this size is probably enough if all we want to do is visualise in igv 
# also probably sufficient res for most coverage comparisons
## (use for IGV)
bamCoverage \
--bam $bam_dir/${ID}.bam \
-o ${bam_dir}_bigWigs_deeptools/${ID}.10.bw \
--binSize 10 \
--normalizeUsing CPM \
--extendReads \
-p 2

# also make 100bp tile coverage plots to compare to other 100bp tile data (eg WGBS)
# output as bed file
# skip zeros (can impute later if needed)
bamCoverage \
--bam $bam_dir/${ID}.bam \
-o ${bam_dir}_100bp_tiles/${ID}.bed \
--binSize 100 \
--normalizeUsing CPM \
--extendReads \
--skipNAs \
--outFileFormat bedgraph \
-p 2

############ get frag sizes ##############
# get frag sizes to check level of tagmentation using picard
# path to call picard jar - eg  - "/home/uqpcrisp/software/picard.jar"

########## Modules #################
module load java/11.0.18
# picard needs R
module load r/4.2.1-foss-2022a

java -jar ~/software/picard.jar CollectInsertSizeMetrics \
I=$bam_dir/${ID}_sorted.bam \
O=${bam_dir}_insert_metrics/${ID}_insert_size_metrics.txt \
H=${bam_dir}_insert_metrics/${ID}_insert_size_histogram.pdf \
M=0.5

echo finished summarising