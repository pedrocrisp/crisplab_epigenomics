#!/usr/bin/Rscript
##########

# 2022/10/27
# Peter Crisp
# Potential batch script to run CSAW outputs to make multi-sample heatmaps 
# User guide: http://bioconductor.org/books/3.15/csawBook/ 

# NOTE: currently hard coded for paired-end sequencing

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
contrast <- args[1]
contrast
region_file <- args[2]
path_to_data_files <- args[3]
path_to_data_files
blacklist <- args[4]
blacklist
filter_FC <- args[5]
filter_FC
bin_size <- args[6]
bin_size
window_spacing <- args[7]
window_spacing

######## de bug
# # args
# 
# # test directors
# qsub -I -l walltime=24:00:00,nodes=1:ppn=4,mem=40gb -A UQ-SCI-SAFS
# cd /scratch/project/crisp002/vanessa/maize_tissues_UMRseq_pete/analysis
# module load R/3.5.0-gnu
# 
# # load R
# R
# # blacklist
# blacklist = "tissues_bams/B73_leaf3_V3_gDNA_test_peaks.narrowPeak"
# # data
# path_to_data_files="tissues_bams"
# 
# # R -f ~/gitrepos/crisplab_epigenomics/methylome/11a-make_DMR_test_table.R --args samples_csaw.txt
# # DMR_contrasts_table_file="DMR_tests_combos_all_table.tsv"
# contrast = "B73_leaf3_V3.vs.Mo17_leaf3_V3"
# 
# filter_FC = 3
# bin_size = 100
# window_spacing = 50
# 
# region_file = "/scratch/project/crisp006/pete/UMRseq_Aim2_run1_mach_II/maize_mach_II/analysis/B_M_bams_2_CSAW/B73_leaf3_V3.vs.Mo17_leaf3_V3_FCF3_bin100_space50/B73_leaf3_V3.vs.Mo17_leaf3_V3_differential_UMRs_CSAW_sig_metadata.bed"
###########################
########
###########################
library(tidyverse)
library(csaw)
library(rtracklayer)
library(GenomicRanges)
library(edgeR)
library(statmod)

library("pheatmap")
library("RColorBrewer")
library(matrixStats)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
                           axis.title=element_text(size=8),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.title=element_text(size=8),
                           legend.text=element_text(size=8))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

# text sizes
text_size_theme_8_labels <- theme(axis.text=element_text(size=8),
                                  axis.title=element_text(size=8),
                                  axis.text.x=element_text(angle = 45, hjust = 1),
                                  legend.title=element_text(size=8),
                                  legend.text=element_text(size=8),
                                  panel.background = element_rect(fill = "transparent") # bg of the panel
                                  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
                                  , panel.grid.major = element_blank() # get rid of major grid
                                  , panel.grid.minor = element_blank() # get rid of minor grid
                                  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
                                  , legend.box.background = element_rect(fill = "transparent")) # get rid of legend panel bg
# magic geom_text conversion ratio
# https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
label_size = 25.4/72 * 8

geom.text.size = 8 / (14/5)

####### set up
projectFolder <- paste0(path_to_data_files, "_CSAW")
dir.create(projectFolder)

outFolder <- paste0(projectFolder, "/heatmaps")
dir.create(outFolder)

# CSAW analysis

################
## import files
################

# import blacklist as GRanges object using rtracklayer
gr_obj =  import(blacklist)
blacklist <- gr_obj

# import regiosn file
region_file =  import(region_file) 

# set up samples
# find file ending with ".bam"
bam.files <- dir("tissues_bams", pattern = "*.bam$", full.names = T, ) 
# get file names
bam.files

# or

bam.files <- bam.files[c(1:4, 7:8, 13:14)]
bam.files

# celltype <- c("leaf3v3", "leaf3v3", "col", "col", "emb", "emb", "endo", "endo", "plumule", "plumule", "root", "root", "scutellum", "scutellum")
celltype <- c("leaf3v3", "leaf3v3", "col", "col", "endo", "endo", "scutellum", "scutellum")
celltype

data.frame(BAM=bam.files, CellType=celltype)

param <- readParam(minq=10, discard=blacklist, pe="both", max.frag=500)


############
# get counts for heatmap and also for quant vs on/off filter
############

# Pull counts for windows in the supplied region file
# them TMM normalise the counts

################
################
# First read in all the data into windows (eg 100bp) and TMM normalise to get normfactors and lib sizes
## this step takes a while
win.data <- windowCounts(bam.files, param=param, width=as.double(bin_size), ext=50, spacing=as.double(window_spacing))
win.data

## potential alternative means to calculate library scaling factors:
#  https://www.biostars.org/p/413626/#414440
# starting with the raw counts per window
## edgeR:: calcNormFactors
NormFactor <- calcNormFactors(object = assay(win.data), method = "TMM")

## raw library size:
LibSize <- colSums(assay(win.data))

## calculate size factors (per 10M):
SizeFactors <- NormFactor * LibSize / 10000000

## Reciprocal, please read section below:   
SizeFactors.Reciprocal <- 1/SizeFactors

################
################
# Now read in the data over the regions of interest
# then normalise

# determine counts over regions
out.ranges.sig.counts <- regionCounts(bam.files, regions=region_file, ext=50, param=param)

out.ranges.sig.counts$size <- getWidths(out.ranges.sig.counts)

# give the counts dataframe row/col names
coords=paste(as.vector(seqnames(out.ranges.sig.counts)),as.vector(start(out.ranges.sig.counts)-1),as.vector(end(out.ranges.sig.counts)),sep="_")
size= as.vector(as.vector(end(out.ranges.sig.counts))-as.vector(start(out.ranges.sig.counts)-1))
# take a look
head(colData(out.ranges.sig.counts))
head(rowData(out.ranges.sig.counts))
head(assays(out.ranges.sig.counts)$counts)
head(size)

samples=bam.files

rownames(out.ranges.sig.counts)=coords
colnames(out.ranges.sig.counts)=samples

out.ranges.sig.counts

# make into a tibble for adding metadata and normalising
counts.tbl <- as_tibble(assays(out.ranges.sig.counts)$counts) %>%
  mutate(coords = coords, size = size) %>%
  mutate(size_factor = size/100) %>%
  select(-size) %>%
  pivot_longer(cols = -c(coords, size_factor), names_to = "BAM", values_to = "raw_counts")

counts.tbl

# make a sample key for normalisation
sample_norm_key <- tibble(BAM = samples, SizeFactors.Reciprocal = SizeFactors.Reciprocal, CellType=celltype)

sample_norm_key

# join and normalise and spread
counts.tbl.norm <- counts.tbl %>% left_join(sample_norm_key, by = "BAM") %>%
  mutate(norm_counts = raw_counts * SizeFactors.Reciprocal/size_factor) %>%
  select(-size_factor) %>%
  pivot_wider(names_from = BAM, values_from = norm_counts, id_cols = coords)

counts.tbl.norm

# write the normalised table in case you want to make more or different heatmaps (this processing took a while!)
write_csv(x = counts.tbl.norm, file = paste0(outFolder, "/", contrast, "_heatmap_norm_counts.csv"), col_names = T)

# join and normalise and spread
counts.tbl.norm_average <- counts.tbl %>% left_join(sample_norm_key, by = "BAM") %>%
  mutate(norm_counts = raw_counts * SizeFactors.Reciprocal/size_factor) %>%
  group_by(CellType, coords) %>%
  summarise(norm_counts_average = mean(norm_counts, na.rm = TRUE)) %>%
  ungroup() %>% 
  select(-size_factor) %>%
  pivot_wider(names_from = CellType, values_from = norm_counts_average)

counts.tbl.norm_average

# write the normalised table in case you want to make more or different heatmaps (this processing took a while!)
write_csv(x = counts.tbl.norm_average, file = paste0(outFolder, "/", contrast, "_heatmap_norm_counts_average.csv"), col_names = T)


################
################
#### heatmap time

# # local debug
# counts.tbl.norm <- read_csv("tissues_bams_CSAW/heatmaps/B73_leaf3_V3.vs.Mo17_leaf3_V3_heatmap_norm_counts.csv")
# counts.tbl.norm_average <- read_csv("tissues_bams_CSAW/heatmaps/B73_leaf3_V3.vs.Mo17_leaf3_V3_heatmap_norm_counts_average.csv")
# contrast = "B73_leaf3_V3.vs.Mo17_leaf3_V3"
# outFolder <- ("tissues_bams_CSAW/heatmaps")

max_rows_to_include = 25000

save_pheatmap_pdf_size <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png_size <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png_size_cm <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height, res = 400, units = "cm")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

##############
# sample counts
counts.tbl.norm %>% pivot_longer(names_to = "tissues", cols = -coords) %>% filter(value != 0) %>% summarise(min = min(value), max = max(value))

counts.tbl.norm.log <- counts.tbl.norm %>% 
  pivot_longer(names_to = "tissues", cols = -coords) %>% 
  mutate(log_value = log2(value)) %>%
  mutate(log_value = ifelse(is.nan(log_value) | is.infinite(log_value), -8, log_value)) %>%
  pivot_wider(names_from = tissues, values_from = log_value, -value)
  
counts.tbl.norm.log

heatmap_input <- counts.tbl.norm.log

# hcluts doesnt linke more than 50,000 rows - so better to subsample
frac_50K <- max_rows_to_include/length(heatmap_input$coords)

big_data_summary_wide_subset <- heatmap_input

if(frac_50K < 1){
  set.seed(27)
  big_data_summary_wide_subset <- big_data_summary_wide_subset %>% sample_frac(frac_50K)
}

big_data_summary_wide_subset

heatmap_data <- data.frame(big_data_summary_wide_subset)
row.names(heatmap_data) <- heatmap_data[,1]
heatmap_data <- data.matrix(heatmap_data[,-c(1)])

# probably put the samples in a better order?

pheat_colour = rev(brewer.pal(11,"RdBu"))

# distance.row = dist(heatmap_data, method = "euclidean")
# cluster.row = hclust(distance.row, method = "ward.D")

h=30
w=12

p <- pheatmap(heatmap_data, silent = T, cluster_rows=TRUE, cluster_cols=T, color = pheat_colour, show_rownames = F,  show_colnames = TRUE, fontsize = 10, na_col = "purple")
save_pheatmap_pdf_size(p, paste0(outFolder, "/", contrast, "_pheatmap", "_max_", max_rows_to_include, "_sampleNormCounts",  ".pdf"), height = h/2.54, width = w/2.54)
save_pheatmap_png_size_cm(p, paste0(outFolder, "/", contrast, "_pheatmap", "_max_", max_rows_to_include, "_sampleNormCounts", ".png"), height = h, width = w)

##############
# average counts
counts.tbl.norm_average %>% pivot_longer(names_to = "tissues", cols = -coords) %>% filter(value != 0) %>% summarise(min = min(value), max = max(value))

counts.tbl.norm.log <- counts.tbl.norm_average %>% 
  pivot_longer(names_to = "tissues", cols = -coords) %>% 
  mutate(log_value = log2(value)) %>%
  mutate(log_value = ifelse(is.nan(log_value) | is.infinite(log_value), -8, log_value)) %>%
  pivot_wider(names_from = tissues, values_from = log_value, -value)

counts.tbl.norm.log

heatmap_input <- counts.tbl.norm.log

# hcluts doesnt linke more than 50,000 rows - so better to subsample
frac_50K <- max_rows_to_include/length(heatmap_input$coords)

big_data_summary_wide_subset <- heatmap_input

if(frac_50K < 1){
  set.seed(27)
  big_data_summary_wide_subset <- big_data_summary_wide_subset %>% sample_frac(frac_50K)
}

big_data_summary_wide_subset

heatmap_data <- data.frame(big_data_summary_wide_subset)
row.names(heatmap_data) <- heatmap_data[,1]
heatmap_data <- data.matrix(heatmap_data[,-c(1)])

# probably put the samples in a better order?

pheat_colour = rev(brewer.pal(11,"RdBu"))

# distance.row = dist(heatmap_data, method = "euclidean")
# cluster.row = hclust(distance.row, method = "ward.D")

h=30
w=12

p <- pheatmap(heatmap_data, silent = T, cluster_rows=TRUE, cluster_cols=T, color = pheat_colour, show_rownames = F,  show_colnames = TRUE, fontsize = 10, na_col = "purple")
save_pheatmap_pdf_size(p, paste0(outFolder, "/", contrast, "_pheatmap", "_max_", max_rows_to_include, "_averageNormCounts",  ".pdf"), height = h/2.54, width = w/2.54)
save_pheatmap_png_size_cm(p, paste0(outFolder, "/", contrast, "_pheatmap", "_max_", max_rows_to_include, "_averageNormCounts", ".png"), height = h, width = w)












###############
print(paste0("Phew, got to the end comparing ", sample1, " vs ", sample2))
###############