#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N merge
#PBS -r y
#PBS -m abej
#PBS -M p.crisp@uq.edu.au

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: array_ID is ${PBS_ARRAY_INDEX}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

########## Modules #################
module load samtools/1.9

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

cd analysis

# created a folder called bsmapped_filtered_merge
mkdir -p bsmapped_filtered_merge

if [ "$paired_end" == "yes" ]
then

########## Run #################

# merge
samtools merge \
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap.bam \
bsmapped_filtered/${ID}_${split_1_name}_sorted_MarkDup_pairs_clipOverlap.bam \
bsmapped_filtered/${ID}_${split_2_name}_sorted_MarkDup_pairs_clipOverlap.bam

# grab header for later
samtools view -H bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap.bam > \
bsmapped_filtered_merge/${ID}_header.sam

# view to sam (can merge output a sam? - yes)
samtools view bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap.bam \
-o bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap.sam

# Add NH flag
# This took around 60Gb and an hour
awk -v FS='\t' -v OFS='\t' 'NR == FNR {T[$1]++; next} {print $0, "NH:i:"T[$1]}' \
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap.sam
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap.sam > \
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH.sam

# Count NH flag and sam flag distribution
grep -oP 'NH:i:\d+' \
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH.sam | \
sed 's/:/\t/g' | cut -f3 | sort | uniq -c > \
bsmapped_filtered_merge/${ID}_merged_NH_distribution.txt

# sam flag distribution
# This takes a couple if gig and 20 min or so
cut -f 2 bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH.sam | \
sort | uniq -c > \
bsmapped_filtered_merge/${ID}_merged_sam_flag_distribution.txt

# Filter to reads (QNAME) occurring twice in sam
# Because these are properly paired reads and QNAME marks reads from the same
# molecule/fragment PE reads have the same QNAME. I have previously filtered to
# mapped and properly paired so filtering noe to NH == 2 should get us what we want...
grep -w 'NH:i:2' \
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH.sam > \
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH_unique.sam

# Add header, sort compress
cat bsmapped_filtered_merge/${ID}_header.sam \
bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH_unique.sam | \
samtools sort \
-O bam \
-@ 2 \
-m 2G \
-o bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH_unique_sorted.bam

# Put his file back in the bsmapped folder
# Use original sample name and giv suffix consistent with the pipeline.
mv -rv bsmapped_filtered_merge/${ID}_merged_sorted_MarkDup_pairs_clipOverlap_NH_unique_sorted.bam \
bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

# once this completes successfully delete all bams/sams in the merge folder (or uncomment this)
# rm -rv bsmapped_filtered_merge/*bam
# rm -rv bsmapped_filtered_merge/*sam

elif [ "$paired_end" == "no" ]
then

echo "uh oh - havent thought about how to deal with SE yet... better write some more code"

else

echo "library type not specified correctly, please indicate yes for PE or no for SE"

fi
