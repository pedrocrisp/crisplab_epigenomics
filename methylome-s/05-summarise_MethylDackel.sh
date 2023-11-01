#!/bin/bash -l
#SBATCH --job-name summarise_MethylDackel
#SBATCH --requeue
#SBATCH --partition=general

# for barley increase walltime to 48hr + 80 Gb

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
conda activate $conda_enviro
module load bedtools/2.30.0-gcc-10.3.0
module load samtools/1.13-gcc-10.3.0

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p MethylDackel
mkdir -p MethylDackel_bigwigs
mkdir -p ConversionRate
mkdir -p MethylDackel_mbias

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
-o MethylDackel/${ID}_methratio_head \
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
bwa-meth-filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam \
MethylDackel_mbias/${ID}_methratio_mbias

# now make a bigwig
#Make bigWigs per context

## remove header
#tail -n+2 MethylDackel/${ID}_methratio_head_CpG.bedGraph > MethylDackel/${ID}_methratio_CG.bedGraph
#tail -n+2 MethylDackel/${ID}_methratio_head_CHG.bedGraph > MethylDackel/${ID}_methratio_CHG.bedGraph
#tail -n+2 MethylDackel/${ID}_methratio_head_CHH.bedGraph > MethylDackel/${ID}_methratio_CHH.bedGraph

# make 6 columns and no header for possible dowstream analysis and for claculating conversion rate (can delete after calculating conversion rate)
# NOTE: this whole genome methratio.bedGraph files only have % meethylated nit the conuts of C and CT per cytosine - this is only kept in the per chromosome files below
# case-sensitive sort
# printf '%s\n' B A b a | LC_COLLATE=en_US.UTF-8 sort
# printf '%s\n' B A b a | LC_COLLATE=C sort
# LC_COLLATE=C

awk -F$"\\t" \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6}' \
MethylDackel/${ID}_methratio_head_CpG.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n - > MethylDackel/${ID}_methratio_CG.bedGraph

awk -F$"\\t" \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6}' \
MethylDackel/${ID}_methratio_head_CHG.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n - > MethylDackel/${ID}_methratio_CHG.bedGraph

awk -F$"\\t" \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6}' \
MethylDackel/${ID}_methratio_head_CHH.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n - > MethylDackel/${ID}_methratio_CHH.bedGraph

# make 4 columns and no header for bedGraphToBigWig
# bw cant pipe to bedgraph to bigwig...
mkdir -p tmp_bedgraphs

awk -F$"\\t" \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4}' \
MethylDackel/${ID}_methratio_head_CpG.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n - > tmp_bedgraphs/${ID}_methratio_CG.bedGraph

awk -F$"\\t" \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4}' \
MethylDackel/${ID}_methratio_head_CHG.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n - > tmp_bedgraphs/${ID}_methratio_CHG.bedGraph

awk -F$"\\t" \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4}' \
MethylDackel/${ID}_methratio_head_CHH.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n - > tmp_bedgraphs/${ID}_methratio_CHH.bedGraph

# make the bedgraphs
bedGraphToBigWig "tmp_bedgraphs/${ID}_methratio_CG.bedGraph" ${chrom_sizes_file} \
"MethylDackel_bigwigs/${ID}_MethylDackel_CG.bigWig"
bedGraphToBigWig "tmp_bedgraphs/${ID}_methratio_CHG.bedGraph" ${chrom_sizes_file} \
"MethylDackel_bigwigs/${ID}_MethylDackel_CHG.bigWig"
bedGraphToBigWig "tmp_bedgraphs/${ID}_methratio_CHH.bedGraph" ${chrom_sizes_file} \
"MethylDackel_bigwigs/${ID}_MethylDackel_CHH.bigWig"

# delete the tmp bedgraphs
rm -rv tmp_bedgraphs

# Now split by chromosome - uncomment if you want to split into per Chr files
# split bedGraph by chromosome
# note - you should probably remove contigs before this step or you could end up with 1000s of files
# eg samtools faidx Vfaba.v1.fasta chr1L chr1S chr2 chr3 chr4 chr5 chr6 > Vfaba.v1_no_contigs.fasta

if [ "$chr_or_genome" == "chromosome" ]
then

awk -F$"\\t" -v ID=$ID \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6 > "MethylDackel/"ID"_"$1"_methratio_CG.bedGraph"}' \
MethylDackel/${ID}_methratio_head_CpG.bedGraph

awk -F$"\\t" -v ID=$ID \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6 > "MethylDackel/"ID"_"$1"_methratio_CHG.bedGraph"}' \
MethylDackel/${ID}_methratio_head_CHG.bedGraph

awk -F$"\\t" -v ID=$ID \
'BEGIN {OFS = FS} (NR>1){print $1, $2, $3, $4, $5, $6 > "MethylDackel/"ID"_"$1"_methratio_CHH.bedGraph"}' \
MethylDackel/${ID}_methratio_head_CHH.bedGraph
## If you wanted to Add context - modify above:
# {print $1, $2, $3, $4, $5, $6, "CG" > "MethylDackel/"ID"_"$1"_methratio_CG.bedGraph"}
# etc

fi

# rm files with header
rm -rv MethylDackel/${ID}_methratio_head*

########################

# I cant remember why i commented this out - maybe when this was in dev, the faba genome didnt have a Chl and the data was published so I didnt bother (and couldnt not QC?)
# i think i should mandate this though, so better to make this work next time you run this pipeline and mandate a chloroplast for conversion eff calc!!!

#ChrC_name=ChrC

## conversion rate
#awk -F$"\\t" \
#'BEGIN {OFS = FS} {sum1 += $6; sum2 +=$5} END {print sum1, sum2 , 100-((sum2/sum1)*100)}' \
#MethylDackel/${ID}_${ChrC_name}_MethylDackel.bedGraph > ConversionRate/${ID}_conversion_rate.txt
## CT  C Conversion_rate
## 365     2       99.4521


echo finished summarising
