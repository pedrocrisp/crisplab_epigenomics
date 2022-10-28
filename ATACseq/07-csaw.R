#!/usr/bin/Rscript
##########

# 2022/10/21
# Peter Crisp
# Batch script to run CSAW on sample pairs
# User guide: http://bioconductor.org/books/3.15/csawBook/ 

# NOTE: currently hard coded for paired-end sequencing

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
contrast <- args[1]
contrast
DMR_contrasts_table_file <- args[2]
DMR_contrasts_table_file
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
# cd /scratch/project/crisp006/pete/UMRseq_Aim2_run1_mach_II/maize_mach_II/analysis
# module load R/3.5.0-gnu
# 
# # load R
# R
# # blacklist
# blacklist = "B_M_bams/B73_leaf3_V3_gDNA_test_peaks.narrowPeak"
# # data
# path_to_data_files="B_M_bams"
# 
# # R -f ~/gitrepos/crisplab_epigenomics/methylome/11a-make_DMR_test_table.R --args samples_csaw.txt
# DMR_contrasts_table_file="DMR_tests_combos_all_table.tsv"
# contrast = "B73_leaf3_V3.vs.Mo17_leaf3_V3"
# 
# filter_FC = 3
# bin_size = 100
# window_spacing = 50
###########################
########
###########################
library(tidyverse)
library(csaw)
library(rtracklayer)
library(GenomicRanges)
library(edgeR)
library(statmod)

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

outFolder <- paste0(projectFolder, "/", contrast, "_FCF", filter_FC, "_bin", bin_size, "_space", window_spacing)
dir.create(outFolder)


DMR_contrasts_table <- read_tsv(DMR_contrasts_table_file)

sample1 = pull(filter(DMR_contrasts_table, test_name == contrast), sample1)
sample1

sample2 = pull(filter(DMR_contrasts_table, test_name == contrast), sample2)
sample2

print(paste0("Calling differential UMRs for ", sample1, " vs ", sample2))

# CSAW analysis

################
## import files
################

# import blacklist as GRanges object using rtracklayer
gr_obj =  import(blacklist)
blacklist <- gr_obj

# set up samples
# find file ending with ".bam"
bam.files <- c(dir(paste0(path_to_data_files, "/"), pattern = paste0("^", sample1, ".*.bam$"), full.names = T, ),
               dir(paste0(path_to_data_files, "/"), pattern = paste0("^", sample2, ".*.bam$"), full.names = T, ))

# get file names
bam.files

# this uses regex - be careful you dont have any file names that com eup in multiple regex... eg "leaf" and "lead_v3" would be a problem

# bam.files <- bam.files[c(1,2,4,5)]
# bam.files

sample1_reps <- length(dir(paste0(path_to_data_files, "/"), pattern = paste0("^", sample1, ".*.bam$"), full.names = T, ))

sample2_reps <- length( dir(paste0(path_to_data_files, "/"), pattern = paste0("^", sample2, ".*.bam$"), full.names = T, ))

celltype <- c(rep(sample1, sample1_reps), rep(sample2, sample2_reps))
celltype

data.frame(BAM=bam.files, CellType=celltype)

# this is the setting for paired end sequencing
param <- readParam(minq=10, discard=blacklist, pe="both", max.frag=500)

# ## examine frag length
# ## This take 5 mins or so - so i have commented it out for now to save time in debugging - UNCOMMENT TO RUN once rest of code is good
# x <- correlateReads(bam.files, param=reform(param, dedup=TRUE)) 
# frag.len <- which.max(x) - 1
# frag.len
# # [1] 9 
# 
# pdf(file = paste0(outFolder, "/frag_length.pdf"))
# plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l") 
# abline(v=frag.len, col="red")
# text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")
# dev.off()

################
## Counts per window and filter
################

# Counting reads into windows
# 100bp windows (just because that what we usually do) (default is 150 for nucleosome length)
# not sure what to use for average frag length because this isnt ChIP - using 50 for now, probably need to tune this variable or just use 0
# default window spacing is 50bp
## this step takes a while
win.data <- windowCounts(bam.files, param=param, width=as.double(bin_size), ext=50, spacing=as.double(window_spacing))
win.data

# Filter windows - these are the setting in the example, this filters windows that have less than 3x the average window coverage - seems reasonable
# may need to tweak

# this step takes a while
bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param) 
# identify windows to filter
filter.stat <- filterWindows(win.data, bins, type="global")
min.fc <- as.double(filter_FC)
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
#    Mode   FALSE    TRUE 
# logical 1421853 2906364 
# we keep more than we discard
# consider raising the filter later

# make plot to visualise the filtering
pdf(file = paste0(outFolder, "/window_filter.pdf"))
hist(filter.stat$back.abundances, main="", breaks=50, xlab="Background abundance (log2-CPM)")
threshold <- filter.stat$abundances[1] - filter.stat$filter[1] + log2(min.fc) 
abline(v=threshold, col="red")
dev.off()

# do the filtering on the RangedSummarizedExperiment
filtered.data <- win.data[keep,]

################
## normalise
################

# Trended biases cannot be removed by scaling methods like TMM normalization (Robinson & Oshlack, 2010), as the amount of scaling required varies with the abundance of the window. Rather, non-linear normalization methods must be used. csaw implements a version of the fast loess method (Ballman et al., 2004) that has been modified to handle count data (Lun & Smyth, 2015). This produces a matrix of offsets that can be used during GLM fitting.
win.ab <- filter.stat$abundances[keep]
head(win.ab)
adjc <- log2(assay(filtered.data)+0.5)
head(adjc)
logfc <- adjc[,1] - adjc[,4]
head(logfc)

# Abundance-dependent trend in the log-fold change between two H3K9ac libraries (mature B over pro-B),  across all windows retained after filtering.
pdf(file = paste0(outFolder, "/Abundance-dependent-trend-bias.pdf"))
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
xlab="Average abundance", ylab="Log-fold change")
dev.off()

# determine offsets for normalisation
offsets <- normOffsets(filtered.data, type="loess", se.out=FALSE) 
head(offsets)

#The effect of non-linear normalization can be visualized with a mean-difference plot comparing the first and last librar- ies. Once the offsets are applied to adjust the log-fold changes, the trend is eliminated from the plot (Figure 4). The cloud of points is also centred at a log-fold change of zero. This indicates that normalization was successful in remov- ing the differences between libraries.
norm.adjc <- adjc - offsets/log(2)
norm.fc <- norm.adjc[,1]-norm.adjc[,4]

pdf(file = paste0(outFolder, "/non-linear-normalization-trended-bias.pdf"))
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5), xlab="Average abundance", ylab="Log-fold change")
dev.off()

# The implicit assumption of non-linear methods is that most windows at each abundance are not DB. Any systematic difference between libraries is attributed to bias and is removed. The assumption of a non-DB majority is reasonable for this data set, given that the cell types being compared are quite closely related. However, it is not appropriate in situ- ations where large-scale DB is expected, as removal of the difference would result in loss of genuine DB. An alternative normalization strategy for these situations will be described later in the CBP analysis.

################
## Statistical modelling of biological variability
################

# experiment design
celltype <- factor(celltype)
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)
design

# Estimating the NB dispersion
y <- asDGEList(filtered.data)
y$offset <- offsets
# this step will take a while
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09609 0.12330 0.16317 0.17495 0.22407 0.27821 

#The NB dispersion trend is visualized in Figure 5 as the biological coefficient of variation (BCV), i.e., the square root of the NB dispersion. Note that only the trended dispersion will be used in the downstream steps – the common and tagwise values are only shown for diagnostic purposes. Specifically, the common BCV provides an overall measure of the variability in the data set, averaged across all windows. Data sets with common BCVs ranging from 10 to 20% are considered to have low variability, i.e., counts are highly reproducible. The tagwise BCVs should also be dispersed above and below the fitted trend, indicating that the fit was successful.

#For most data sets, one would expect to see a trend that decreases to a plateau with increasing average abundance. This reflects the greater reliability of large counts, where the effects of stochasticity and technical artifacts (e.g., mapping errors, PCR duplicates) are averaged out. In Figure 5, the range of abundances after filtering is such that the plateau has already been reached. This is still a satisfactory result, as it indicates that the retained windows have low variability and more power to detect DB.

png(file = paste0(outFolder, "/Abundance-dependent-trend-in-the-BCV-for-each-window.png"))
plotBCV(y)
dev.off()

# Estimating the QL dispersion.
# Additional modelling is provided with the QL methods in edgeR (Lund et al., 2012). This introduces a QL dispersion parameter for each window, which captures variability in the NB dispersion around the fitted trend for each window. Thus, the QL dispersion can model window-specific variability, whereas the NB dis- persion trend is averaged across many windows. However, with limited replicates, there is not enough information for each window to stably estimate the QL dispersion. This is overcome by sharing information between windows with empirical Bayes (EB) shrinkage. The instability of the QL dispersion estimates is reduced by squeezing the estimates towards an abundance-dependent trend (Figure 6).

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 13.18   13.33   13.33   13.31   13.33   13.33 

# why are these numbers all 13.33??????

png(file = paste0(outFolder, "/Effect-of-EB-shrinkage-on-the-raw-QL-dispersion-estimate-for-each-window.png"))
plotQLDisp(fit)
dev.off()

# Examining the data with MDS plots. Multi-dimensional scaling (MDS) plots can be used to examine the similarities between libraries. The distance between a pair of libraries on this plot represents the overall log-fold change between those libraries. Ideally, replicates should cluster together while samples from different conditions should be separate. In Figure 7, strong separation in the first dimension is observed between libraries from different cell types. This indicates that significant differences are likely to be present between cell types in this data set.

pdf(file = paste0(outFolder, "/mds.pdf"))
plotMDS(norm.adjc, labels=celltype, col=c("red", "blue")[as.integer(celltype)])
dev.off()

################
## Testing for differenital UMRs ("DB") and controlling the FDR
################

# Testing for DB with QL F-tests. Each window is tested for significant differences between cell types using the QL F-test (Lund et al., 2012). This is superior to the likelihood ratio test that is typically used for GLMs, as the QL F-test accounts for the uncertainity in dispersion estimation. One p-value is produced for each window, representing the evidence against the null hypothesis (i.e., that no DB is present in the window). For this analysis, the comparison is parametrized such that the reported log-fold change for each window represents that of the coverage in pro-B cells over their mature B counterparts.

# to parse the sample name to makeContrasts
# euqivalent to: contrast <- makeContrasts(sample1-sample2, levels=design)
mycontrast = paste0(sample1, "-", sample2) 
cmd <- paste("tmp <- makeContrasts(", mycontrast, ", levels = design)", sep = '"')
contrastEval <- eval(parse(text = cmd))

res <- glmQLFTest(fit, contrast=contrastEval)
head(res$table)

#       logFC     logCPM          F    PValue
# 1 1.1196548 -0.7306550 1.65215581 0.2177414
# 2 0.8588430 -0.5071920 1.23349501 0.2838417
# 3 1.0672588 -0.4065665 2.02976546 0.1742798
# 4 1.4593595 -0.8647425 1.79073457 0.2005337
# 5 0.3156714 -0.9548798 0.09678204 0.7599154
# 6 0.3707627 -0.8507514 0.14431459 0.7092391

# region level FDR correction
# Controlling the FDR across regions. One might attempt to control the FDR by applying the Benjamini-Hochberg (BH) method to the window-level p-values (Benjamini & Hochberg, 1995). However, the features of interest are not windows, but the genomic regions that they represent. Control of the FDR across windows does not guarantee control of the FDR across regions (Lun & Smyth, 2014). The latter is arguably more relevant for the final interpretation of the results.
# Control of the region-level FDR can be provided by aggregating windows into regions and combining the p-values. Here, adjacent windows less than 100 bp apart are aggregated into clusters. Each cluster represents a genomic region. Smaller values of tol allow distinct marking events to kept separate, while larger values provide a broader perspective, e.g., by considering adjacent co-regulated sites as a single entity. Chaining effects are mitigated by setting a maximum cluster width of 5 kbp.
merged <- mergeWindows(rowRanges(filtered.data),tol=100, max.width=5000)

# A combined p-value is computed for each cluster using the method of Simes (1986), based on the p-values of the con- stituent windows. This represents the evidence against the global null hypothesis for each cluster, i.e., that no DB exists in any of its windows. Rejection of this global null indicates that the cluster (and the region that it represents) contains DB. Applying the BH method to the combined p-values allows the region-level FDR to be controlled.
tabcom <- combineTests(merged$id, res$table)
head(tabcom)
# DataFrame with 6 rows and 6 columns
#    nWindows  logFC.up logFC.down             PValue               FDR
#   <integer> <integer>  <integer>          <numeric>         <numeric>
# 1         4         4          0  0.283841742070597 0.618275979267929
# 2         2         0          0  0.759915396348918 0.994300332636627
# 3         1         1          0  0.272135330351028 0.603019010588037
# 4        18        13          0  0.362269729381117 0.710508213518741
# 5        13        13          0 0.0559986698401485 0.208810012625004
# 6         1         0          0  0.909343338491651 0.999970818301873
#     direction
#   <character>
# 1          up
# 2          up
# 3          up
# 4          up
# 5          up
# 6          up

# Each row of the output table contains the statistics for a single cluster, including the combined p-value before and after the BH correction. The nWindows field describes the total number of windows in the cluster. The logFC.up and logFC.down fields describe the number of windows with a log-fold change above 0.5 or below -0.5 in each cluster, respectively. This can be used to determine the direction of DB in each cluster.

# Examining the scope and direction of DB. The total number of DB regions at a FDR of 5% can be easily calculated. 
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
#    Mode   FALSE    TRUE 
# logical  330236   65446 

# Determining the direction of DB is more complicated, as clusters could potentially contain windows that are changing in opposite directions. One approach is to define the direction based on the number of windows changing in each direction, as described above. Another approach is to use the log-fold change of the most significant window as a proxy for the log-fold change of the cluster. This is generally satisfactory, though it will not capture multiple changes in opposite directions. It also tends to overstate the change in each cluster.
tabbest <- getBestTest(merged$id, res$table)
head(tabbest)
# DataFrame with 6 rows and 6 columns
#        best             logFC             logCPM                  F
#   <integer>         <numeric>          <numeric>          <numeric>
# 1         3  1.06725879366625 -0.406566496977074   2.02976546189978
# 2         6 0.370762697956262 -0.850751352993566  0.144314587458563
# 3         7  1.06236485558409  -0.87973467429518   1.29760474569143
# 4         9  1.36889231792817  0.547138234131017   5.66823568999084
# 5        36  1.62114096693505  0.072084369865961   6.95210139352563
# 6        39  0.11665128776015   -1.0188905866489 0.0134009348400389
#              PValue               FDR
#           <numeric>         <numeric>
# 1 0.697119131797429                 1
# 2                 1                 1
# 3 0.272135330351028 0.672619928814318
# 4 0.551576575031542                 1
# 5 0.239537901652749 0.619550088776878
# 6 0.909343338491651                 1

# In the above table, each row contains the statistics for each cluster. Of interest are the best and logFC fields. The former is the index of the window that is the most significant in each cluster, while the latter is the log-fold change of that window. This can be used to obtain a summary of the direction of DB across all clusters/regions.

is.sig.pos <- (tabbest$logFC > 0)[is.sig]
summary(is.sig.pos)
#    Mode   FALSE    TRUE 
# logical   41524   23922 


################
## saving results
################

# Results can be saved to file prior to further manipulation. One approach is to store all statistics in the metadata of a GRanges object. This is useful as it keeps the statistics and coordinates together for each cluster, avoiding problems with synchronization in downstream steps. The midpoint and log-fold change of the best window are also stored.
out.ranges <- merged$region

elementMetadata(out.ranges) <- data.frame(tabcom,
   best.pos=mid(ranges(rowRanges(filtered.data[tabbest$best]))), 
   best.logFC=tabbest$logFC) 

# # uncomment to save the RDS file
#########
# saveRDS(file=paste0(outFolder, "/", contrast, "_differential_UMRs.rds"), out.ranges)
#########

# For input into other programs like genome browsers, results can be saved in a more conventional format. Here, coordinates of DB regions are saved in BED format via rtracklayer, using a log-transformed FDR as the score.
# simplified <- out.ranges[is.sig]
# simplified$score <- -10*log10(simplified$FDR) 
# export(con=paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW.bed"), object=simplified)

# Saving the RangedSummarizedExperiment objects is also recommended. This avoids the need to re-run the time- consuming read counting steps if parts of the analysis need to be repeated. Similarly, the DGEList object is saved so that the edgeR statistics can be easily recovered.
# # uncomment to save the 
#########
# save(file=paste0(outFolder, "/", contrast, "_objects.Rda"), win.data, bins, y)
#########

### code to make a better bed file
out.ranges

# remove organelles
out.ranges.sig <- out.ranges[!out.ranges@seqnames  %in%  c("Pt", "Mt")]
out.ranges.sig

# filter for significance
out.ranges.sig <- out.ranges.sig[out.ranges.sig$FDR < 0.05]

# note some regions are "mixed"
no_call <- out.ranges.sig[!(out.ranges.sig$direction %in% c("up", "down"))]
no_call

# add direction
out.ranges.sig$name_tmp <- out.ranges.sig$direction
out.ranges.sig$name <- ifelse(out.ranges.sig$direction =="up", sample1, ifelse(out.ranges.sig$direction == "down", sample2, "mixed"))

# size cats
out.ranges.sig$size <- width(out.ranges.sig)
out.ranges.sig$size_cat <- ifelse(out.ranges.sig$size < 300, "small", ifelse(out.ranges.sig$size >=300 & out.ranges.sig$size <900, "med", "large"))

# add the fold change
out.ranges.sig$score <- out.ranges.sig$best.logFC

# export better bedfile
export(con=paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_sig_metadata.bed"), object=out.ranges.sig)

# with more metadata
out.ranges.sig.df <- data.frame(out.ranges.sig)
head(out.ranges.sig.df)
write.table( x = out.ranges.sig.df, file = paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_sig_metadata_extend.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

# just get UMRs (up in) in sample1 eg root specific UMRs
out.ranges.sig_sample1 <- out.ranges.sig[out.ranges.sig$name == sample1]
# export better bedfile
export(con=paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_sig_metadata_", sample1, ".bed"), object=out.ranges.sig_sample1)

# just get UMRs (up in) in sample1 eg root specific UMRs
out.ranges.sig_sample2 <- out.ranges.sig[out.ranges.sig$name == sample2]
# export better bedfile
export(con=paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_sig_metadata_", sample2, ".bed"), object=out.ranges.sig_sample2)

out.ranges.sig.small <- out.ranges.sig[out.ranges.sig$size < 300]
export(con=paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_sig_metadata_less_300bp.bed"), object=out.ranges.sig.small)

out.ranges.sig.large <- out.ranges.sig[out.ranges.sig$size >= 300]
export(con=paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_sig_metadata_greater_300bp.bed"), object=out.ranges.sig.large)

############

## analysis of DE UMRs - CSAW
# size distribution
# small, med, large
# histogram

# make into a tibble for tidyverse
out.ranges.sig.tbl <- as_tibble(out.ranges.sig)

out.ranges.sig.tbl

# make summary
size_summary <- out.ranges.sig.tbl %>% mutate(size_cat = factor(size_cat)) %>% group_by(size_cat) %>% summarise(n = n())
size_summary

write_csv(x = size_summary, file = paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_size_summary.csv"), col_names = T)

# all regions
# # A tibble: 3 × 2
#   size_cat      n
#   <fct>     <int>
# 1 large     44556
# 2 med      173413
# 3 small    177523

# significant regions
# # # A tibble: 3 × 2
#   size_cat     n
#   <fct>    <int>
# 1 large    10973
# 2 med      27405
# 3 small    27067

g <- out.ranges.sig.tbl %>%
ggplot(., aes((size))) +
  geom_density() +
  theme_minimal() +
  text_size_theme_8
g
ggsave(plot = g, filename = paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_sig_sizes_histogram.pdf"), h = 4, w = 6)

### logFC stats

FC_summary <- out.ranges.sig.tbl %>%
  mutate(FC_cat = ifelse(best.logFC < 2, "<2", ifelse(best.logFC >=2 & best.logFC <4, "<4", ">=4"))) %>%
  group_by(FC_cat) %>% summarise(n = n())

FC_summary

write_csv(x = FC_summary, file = paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_FC_summary.csv"), col_names = T)

############
# get counts for heatmap and also for quant vs on/off filter
############
################
# Pull counts for significant windows
# 2 possible approaches - https://support.bioconductor.org/p/98672/ 
## 1. extract the counts for the most significant window within the region using getBestTest and the $best
## 2. re-count using regionCounts (this isnt exactly the same but should be  very close)
## Going for option 2

# out.ranges.sig
# win.data <- windowCounts(bam.files, param=param, width=as.double(bin_size), ext=50, spacing=as.double(window_spacing))
# win.data

out.ranges.sig.counts <- regionCounts(bam.files, regions=out.ranges.sig, ext=50, param=param)
coords=paste(as.vector(seqnames(out.ranges.sig.counts)),as.vector(start(out.ranges.sig.counts)-1),as.vector(end(out.ranges.sig.counts)),sep="_")

head(colData(out.ranges.sig.counts))
head(rowData(out.ranges.sig.counts))
head(assays(out.ranges.sig.counts)$counts)

samples=bam.files

rownames(out.ranges.sig.counts)=coords
colnames(out.ranges.sig.counts)=samples

out.ranges.sig.counts

# cpm(assay(out.ranges.sig.counts), lib.size=exp(getOffset(y)))

counts.tbl <- as_tibble(assays(out.ranges.sig.counts)$counts) %>%
  mutate(coords = coords) %>%
  pivot_longer(cols = -coords, names_to = "BAM", values_to = "raw_counts")

counts.tbl

## potential alternative means to calculate library scaling factors:
#  https://www.biostars.org/p/413626/#414440
# starting with the raw counts per window
## edgeR:: calcNormFactors
NormFactor <- calcNormFactors(object = assay(win.data), method = "TMM")

## if you prefer to use the DESeq2 strategy use method="RLE" instead

## raw library size:
LibSize <- colSums(assay(win.data))

## calculate size factors:
SizeFactors <- NormFactor * LibSize / 1000000

## Reciprocal, please read section below:   
SizeFactors.Reciprocal <- 1/SizeFactors

sample_norm_key <- tibble(BAM = samples, SizeFactors.Reciprocal = SizeFactors.Reciprocal, CellType=celltype)

sample_norm_key

# join and normalise and spread
counts.tbl.norm <- counts.tbl %>% left_join(sample_norm_key, by = "BAM") %>%
  mutate(norm_counts = raw_counts * SizeFactors.Reciprocal) %>%
  group_by(CellType, coords) %>%
  summarise(norm_counts_average = mean(norm_counts)) %>%
  ungroup() %>% pivot_wider(names_from = CellType, values_from = norm_counts_average)

counts.tbl.norm

write_csv(x = counts.tbl.norm, file = paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_norm_counts_average.csv"), col_names = T)

# counts.tbl.norm <- read_csv("~/UMRseq-06/analysis/UMRseq_Aim2/UMRseq_Aim2_run1_mach_II/maize_mach_II/analysis/B_M_bams_CSAW/B73_leaf3_V3.vs.Mo17_leaf3_V3_FCF3_bin100_space50/B73_leaf3_V3.vs.Mo17_leaf3_V3_differential_UMRs_CSAW_norm_counts_average.csv")
# sample1 = "B73_leaf3_V3"
# sample2 = "Mo17_leaf3_V3"
# contrast <- "B73vMo17"
# outFolder <- (paste0("csaw/", contrast))

g <- ggplot(counts.tbl.norm, aes(x=log2(B73_leaf3_V3), y=log2(Mo17_leaf3_V3))) + 
  geom_point(size=0.2) +
  theme_minimal() +
  text_size_theme_8
# g
ggsave(plot = g, filename = paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_norm_counts_scatter.pdf"), h = 8, w = 8)
ggsave(plot = g, filename = paste0(outFolder, "/", contrast, "_differential_UMRs_CSAW_norm_counts_scatter2.pdf"), h = 2, w = 2)

# counts.df.norm <- df2$SizeFactors.Reciprocal[match(names(counts.df), df2$sample)][col(counts.df)] * counts.df
# head(counts.df.norm)

###############
print(paste0("Phew, got to the end comparing ", sample1, " vs ", sample2))
###############
