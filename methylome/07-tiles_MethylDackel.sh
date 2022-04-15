#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N tiles_MethylDackel
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

module load R/3.5.0-gnu
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p tiles

########## Run #################

#Run R moudle to:
# 1. sumamrise into 100bp tiles
# 2. fix chr ends in BED file
# 3. output one-based coordinate .txt file for analysis in R
# folder with input files hardcoded "MethylDackel"

for CHROMOSOME in (cat $chromosome_list); do
R -f ~/gitrepos/crisplab_epigenomics/methylome/07-tiles_MethylDackel_tile_mod.R \
--args ${ID}_${CHROMOSOME} MethylDackel ${reference_tile_file_folder} $CHROMOSOME
done

#make bedGraph by sorting and removing cols 4 and 5 with awk

# make bedGraph by removing cols 4 and 5 and sorting
# cut -f1-3,6-6 ./tiles/${ID}_BSMAP_out.txt.100.CG.fixed.bed | sort -k1,1 -k2,2n > ./tiles/${ID}_BSMAP_out.txt.100.CG.fixed.sorted.bg
# cut -f1-3,6-6 ./tiles/${ID}_BSMAP_out.txt.100.CHG.fixed.bed | sort -k1,1 -k2,2n > ./tiles/${ID}_BSMAP_out.txt.100.CHG.fixed.sorted.bg
# cut -f1-3,6-6 ./tiles/${ID}_BSMAP_out.txt.100.CHH.fixed.bed | sort -k1,1 -k2,2n > ./tiles/${ID}_BSMAP_out.txt.100.CHH.fixed.sorted.bg

#Make bigWigs
# At some point soft code the reference, make variable in script call
# To mkae the chrom_sizes file use seqtk compressseqtk comp Zea_mays.AGPv4.dna.toplevel.fa
# seqtk comp Zea_mays.AGPv4.dna.toplevel.fa | cut -f 1-2 > Zea_mays.AGPv4.dna.toplevel.chrom.sizes

# note - on some genomes this bigwig fails because it is not natural sorted...
# but i generally do use these files so no worries...
# maybe somepoint add "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again

# bedGraphToBigWig "./tiles/${ID}_BSMAP_out.txt.100.CG.fixed.sorted.bg" ${chrom_sizes} \
# "./tiles/${ID}_BSMAP_out.txt.100.CG.bigWig"
# bedGraphToBigWig "./tiles/${ID}_BSMAP_out.txt.100.CHG.fixed.sorted.bg" ${chrom_sizes} \
# "./tiles/${ID}_BSMAP_out.txt.100.CHG.bigWig"
# bedGraphToBigWig "./tiles/${ID}_BSMAP_out.txt.100.CHH.fixed.sorted.bg" ${chrom_sizes} \
# "./tiles/${ID}_BSMAP_out.txt.100.CHH.bigWig"

# remove .bg - not required
# meh purge it in step 09 instead
# rm -rv tiles/${ID}*.bg

echo finished summarising
