#!/bin/bash -l
#SBATCH --job-name coverage
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
conda activate $conda_enviro
# dependencies
# conda env with deeptools-hacked and mosdepth installed
# conda create deeptools-hacked python=3.7
# conda activate deeptools-hacked
# conda install -c bioconda mosdepth
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
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p $out_dir

########## Run #################

# deeptools to make coverage matrix
# summarise 2kb either side of feature
# -p 8 threds
# dont skip regions with only zeros; removed: --skipZeros
# treat missing data as zero

########### CG  ############
computeMatrix scale-regions \
-R ${region_dir}/${ID}${suffix} \
-S $CG_bigWig \
-b 2000 -a 2000 \
-p 2 \
-bs 100 \
--regionBodyLength 1000 \
--sortRegions descend \
-o ${out_dir}/${ID}_CG.mat.gz \
--outFileNameMatrix ${out_dir}/${ID}_CG.mat_values.tab \
--outFileSortedRegions ${out_dir}/${ID}_CG.region_order.bed

#plot
plotProfile \
-m ${out_dir}/${ID}_CG.mat.gz \
-out ${out_dir}/${ID}_CG.mat_metaplot.png \
--numPlotsPerRow 1 \
--yMin 0 \
--yMax 100

# heatmap
plotHeatmap \
-m ${out_dir}/${ID}_CG.mat.gz \
-out ${out_dir}/${ID}_CG.mat_heatmap.png \
--interpolationMethod nearest \
--dpi 1200 \
--zMin 0 \
--zMax 100 \
--colorMap RdYlBu
#--colorMap RdYlBu_r

# heatmap
plotHeatmap \
-m ${out_dir}/${ID}_CG.mat.gz \
-out ${out_dir}/${ID}_CG.mat_heatmap_bilinear.png \
--interpolationMethod bilinear \
--dpi 1200 \
--zMin 0 \
--zMax 100 \
--colorMap RdYlBu


########### CHG  ############
computeMatrix scale-regions \
-R ${region_dir}/${ID}${suffix} \
-S $CHG_bigWig \
-b 2000 -a 2000 \
-p 2 \
-bs 100 \
--regionBodyLength 1000 \
--sortRegions descend \
-o ${out_dir}/${ID}_CHG.mat.gz \
--outFileNameMatrix ${out_dir}/${ID}_CHG.mat_values.tab \
--outFileSortedRegions ${out_dir}/${ID}_CHG.region_order.bed

#plot
plotProfile \
-m ${out_dir}/${ID}_CHG.mat.gz \
-out ${out_dir}/${ID}_CHG.mat_metaplot.png \
--numPlotsPerRow 1 \
--yMin 0 \
--yMax 100

# heatmap
plotHeatmap \
-m ${out_dir}/${ID}_CHG.mat.gz \
-out ${out_dir}/${ID}_CHG.mat_heatmap.png \
--interpolationMethod nearest \
--dpi 1200 \
--zMin 0 \
--zMax 100 \
--colorMap RdYlBu
#--colorMap RdYlBu_r

# heatmap
plotHeatmap \
-m ${out_dir}/${ID}_CHG.mat.gz \
-out ${out_dir}/${ID}_CHG.mat_heatmap_bilinear.png \
--interpolationMethod bilinear \
--dpi 1200 \
--zMin 0 \
--zMax 100 \
--colorMap RdYlBu

########### CHH  ############
computeMatrix scale-regions \
-R ${region_dir}/${ID}${suffix} \
-S $CHH_bigWig \
-b 2000 -a 2000 \
-p 2 \
-bs 100 \
--regionBodyLength 1000 \
--sortRegions descend \
-o ${out_dir}/${ID}_CHH.mat.gz \
--outFileNameMatrix ${out_dir}/${ID}_CHH.mat_values.tab \
--outFileSortedRegions ${out_dir}/${ID}_CHH.region_order.bed

#plot
plotProfile \
-m ${out_dir}/${ID}_CHH.mat.gz \
-out ${out_dir}/${ID}_CHH.mat_metaplot.png \
--numPlotsPerRow 1 \
--yMin 0 \
--yMax 100

# heatmap
plotHeatmap \
-m ${out_dir}/${ID}_CHH.mat.gz \
-out ${out_dir}/${ID}_CHH.mat_heatmap.png \
--interpolationMethod nearest \
--dpi 1200 \
--zMin 0 \
--zMax 100 \
--colorMap RdYlBu
#--colorMap RdYlBu_r

# heatmap
plotHeatmap \
-m ${out_dir}/${ID}_CHH.mat.gz \
-out ${out_dir}/${ID}_CHH.mat_heatmap_bilinear.png \
--interpolationMethod bilinear \
--dpi 1200 \
--zMin 0 \
--zMax 100 \
--colorMap RdYlBu


# Run R module to creat 100bp tile bed file
module load r/4.2.1-foss-2022a
R -f ~/gitrepos/crisplab_epigenomics/methylome-s/09-WGBS-cov_feature_summarise.R \
--args $ID $out_dir

echo finished summarising
