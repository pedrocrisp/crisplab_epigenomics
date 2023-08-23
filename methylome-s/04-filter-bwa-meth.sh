#!/bin/bash -l
#SBATCH --job-name bwa-filter
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

######### Modules #################
# requires java 1.8
conda activate $conda_enviro
module load java/11.0.16
module load bamtools/2.5.2-gcc-10.3.0
module load samtools/1.13-gcc-10.3.0

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

cd analysis
mkdir -p bwa-meth-filtered

if [ "$paired_end" == "yes" ]
then

########## Run #################

        # fix improperly paird reads - specifically discordant read pairs that are incorrectly mraked as concordant by bsmap
        # picard FixMateInformation
        # didnt work...
        #java -jar /home/springer/pcrisp/software/picard.jar FixMateInformation \
        #I=bsmapped/${ID}.bam \
        #O=bwa-meth-filtered/${ID}_sorted.bam

        #remove PCR duplicates, must be sorted by coordinate using pickard

        # bams already co-ordinate sorted by samtools, this step seems unnecessary
        # Also causesing issues: bsmap is reporting PE reads as properly mapped where they hit different chromosomes, solution: skip step
        # uncomment to re-instate

        #java -jar /home/springer/pcrisp/software/picard.jar SortSam \
        #INPUT=bsmapped/${ID}.bam \
        #OUTPUT=bsmapped/${ID}_sorted.bam \
        #SORT_ORDER=coordinate

        #mark duplicates
        #requires sorted input - using samtools sort in bsmap step (co-ordinate sorted)
        # if co-ordinate sorted then pairs where the mate is unmapped or has secondary alignment are not marked as duplicate
        # ASSUME_SORTED=true because sorting performed with samtools but samtools doesnt seem to add this flag to the headder
        java -jar ~/software/picard.jar MarkDuplicates \
        I=bwa-meth/${ID}_sorted.bam \
        O=bwa-meth-filtered/${ID}_sorted_MarkDup.bam \
        METRICS_FILE=bwa-meth-filtered/${ID}_MarkDupMetrics.txt \
        ASSUME_SORTED=true \
        CREATE_INDEX=False \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=true

        # keep properly paired reads using bamtools package
        # note that some reads marked as properly paired by bsmap actually map to different chromosomes
        bamtools filter \
        -isMapped true \
        -isPaired true \
        -isProperPair true \
        -in bwa-meth-filtered/${ID}_sorted_MarkDup.bam \
        -out bwa-meth-filtered/${ID}_sorted_MarkDup_pairs.bam

        # clip overlapping reads using bamUtils package
        bam clipOverlap \
        --in bwa-meth-filtered/${ID}_sorted_MarkDup_pairs.bam \
        --out bwa-meth-filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam \
        --stats

        #index bam
        # index
        samtools index -c -@ 20 bwa-meth-filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

elif [ "$paired_end" == "no" ]
then

# remove duplicate
java -jar ~/software/picard.jar MarkDuplicates \
I=bwa-meth/${ID}_sorted.bam \
O=bwa-meth-filtered/${ID}_sorted_MarkDup.bam \
METRICS_FILE=bwa-meth-filtered/${ID}_MarkDupMetrics.txt \
ASSUME_SORTED=true \
CREATE_INDEX=False \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=true

# rename to keep script naming convention consistent
mv bwa-meth-filtered/${ID}_sorted_MarkDup.bam bwa-meth-filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

# index
samtools index -c -@ 20 bwa-meth-filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam


else

echo "library type not specified correctly, please indicate yes for PE or no for SE"

fi
