#!/bin/bash -l
#SBATCH --job-name tiles_MethylDackel
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

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p tiles

########## Run #################

# Run R module to:
# 1. sumamrise into 100bp tiles
# 2. fix chr ends in BED file
# 3. output one-based coordinate .txt file for analysis in R
# folder with input files hardcoded "MethylDackel"

if [ "$chr_or_genome" == "chromosome" ]
then

for CHROMOSOME in $(cat $chromosome_list); do
R -f ~/gitrepos/crisplab_epigenomics/methylome-s/07-tiles_MethylDackel_tile_mod.R \
--args ${ID}_${CHROMOSOME} MethylDackel ${reference_tile_file_folder_or_file}/${CHROMOSOME}_100bp_tiles_zBased_sites_counts.txt
done

elif [ "$chr_or_genome" == "genome" ]
then

R -f ~/gitrepos/crisplab_epigenomics/methylome-s/07-tiles_MethylDackel_tile_mod.R \
--args ${ID} MethylDackel ${reference_tile_file_folder_or_file}

else

echo "whole or split genome not specificed, please indicate chromosome or genome"

fi

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
