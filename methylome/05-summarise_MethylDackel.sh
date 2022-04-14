#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N summarise_MethylDackel
#PBS -r y
#PBS -m abej
#PBS -M p.crisp@uq.edu.au

# for barley increase walltime to 48hr + 80 Gb

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
conda activate py3.7
module load bedtools
module load samtools/1.9

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p MethylDackel
mkdir -p ConversionRate

########## Run #################
# MethylDackel
# extract CHG and CHH too
# filter sites that have 10% reverse stard 'G' at T's could be mutations, so ignore these if 3 or more reads
# if there is M bias - adjust --nOT 0,0,0,0 and --nOB 0,0,0,0

MethylDackel extract \
--CHG \
--CHH \
-@ $cores \
--minOppositeDepth=3 \
--maxVariantFrac=0.1 \
--nOT 0,0,0,0 \
--nOB 0,0,0,0 \
-o MethylDackel/${ID}_methratio \
${genome_reference} \
bwa-meth-filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

# out put

#    The chromosome/contig/scaffold name
#    The start coordinate
#    The end coordinate
#    The methylation percentage rounded to an integer
#    The number of alignments/pairs reporting methylated bases
#    The number of alignments/pairs reporting unmethylated bases

#track type="bedGraph" description="MethylDackel/SRR5724529_8k_methratio CHG methylation levels"
#Chr1    24540   24541   0       0       1
#Chr1    24547   24548   0       0       1
#Chr1    309269  309270  0       0       1
#Chr1    309280  309281  0       0       1
#Chr1    309332  309333  0       0       1

# make a combined files - dont think we need this... nor does it work?
#sort -k2,2n <(tail -n+2 -q MethylDackel/${ID}_methratio*bedGraph) | awk '{print "Caudatum-SC237-14E\t"$0}' > bwa-meth_filtered/${ID}.allC
#sort -k2,2n <(tail -n+2 -q *bedGraph) | awk '{print "Caudatum-SC237-14E\t"$0}' > Caudatum-SC237-14E.allC

# get bias metrics
MethylDackel mbias \
${genome_reference} \
bwa-meth-filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
MethylDackel/${ID}_methratio_mbias

# now make a bigwig
#Make bigWigs per context
bedGraphToBigWig "MethylDackel/${ID}_methratio_CpG.bedGraph" ${chrom_sizes_file} \
"MethylDackel/${ID}_MethylDackel_CG.bigWig"
bedGraphToBigWig "MethylDackel/${ID}_methratio_CHG.bedGraph" ${chrom_sizes_file} \
"MethylDackel/${ID}_MethylDackel_CHG.bigWig"
bedGraphToBigWig "MethylDackel/${ID}_methratio_CHH.bedGraph" ${chrom_sizes_file} \
"MethylDackel/${ID}_MethylDackel_CHH.bigWig"

## Add context - uncomment if you want to add column with context
#awk -F$"\\t" \
#'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6, "CG"}' \
# MethylDackel/${ID}_methratio_CpG.bedGraph > MethylDackel/${ID}_methratio_CpG_context.bedGraph
#
# awk -F$"\\t" \
# 'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6, "CHG"}' \
#  MethylDackel/${ID}_methratio_CHG.bedGraph > MethylDackel/${ID}_methratio_CHG_context.bedGraph
#
#  awk -F$"\\t" \
#  'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6, "CHH"}' \
#   MethylDackel/${ID}_methratio_CHH.bedGraph > MethylDackel/${ID}_methratio_CHH_context.bedGraph

#

## Now split by chromosome - uncomment if you want to split into per Chr files
## split bedGraph by chromosome
#awk_make_bedGraph_context='BEGIN {OFS = FS} (NR>1){
#  print $1, $2, $3, $4, $5, $6 > "MethylDackel/"ID"_"$1"_MethylDackel.bedGraph"
#}
#'
#
#awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" "MethylDackel/${ID}_methratio_CpG.bedGraph"
#awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" "MethylDackel/${ID}_methratio_CHG.bedGraph"
#awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" "MethylDackel/${ID}_methratio_CHH.bedGraph"

########################

#ChrC_name=ChrC

## conversion rate
#awk -F$"\\t" \
#'BEGIN {OFS = FS} {sum1 += $6; sum2 +=$5} END {print sum1, sum2 , 100-((sum2/sum1)*100)}' \
#MethylDackel/${ID}_${ChrC_name}_MethylDackel.bedGraph > ConversionRate/${ID}_conversion_rate.txt
## CT  C Conversion_rate
## 365     2       99.4521


echo finished summarising
