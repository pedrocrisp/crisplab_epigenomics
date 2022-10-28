#!/usr/bin/Rscript
##########

# 2022/10/28
# Peter Crisp
# Batch script to annotate CSAW output gene annotations
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
# qsub -I -l walltime=24:00:00,nodes=1:ppn=4,mem=40gb -A UQ-SCI-SAFS
# cd /scratch/project/crisp006/pete/UMRseq_Aim2_run1_mach_II/maize_mach_II/analysis
# 
# sample_prefix=B73_leaf3_V3.vs.Mo17_leaf3_V3
# CSAW_settings=FCF3_bin100_space50
# DE_folder=/scratch/project/crisp006/pete/UMRseq_Aim2_run1_mach_II/maize_mach_II/analysis/B_M_bams_2_CSAW
# DE_UMR_bed=${DE_folder}/${sample_prefix}_${CSAW_settings}/${sample_prefix}_differential_UMRs_CSAW_sig_metadata.bed
# 
# outFolder=${DE_folder}/${sample_prefix}_${CSAW_settings}/annotation
# 
# mkdir $outFolder
# 
# echo $DE_UMR_bed
# 
# module load bedtools
# 
# bedtools closest \
# -a $DE_UMR_bed \
# -b ~/refseqs/maize/Zea_mays_AGPv4_36_fixed_introns_gene_ncRNA_synteny_miRbase_space_V_stranded_bed6.bed \
# -mdb all \
# -t all \
# -D b \
# -g ~/refseqs/maize/Zea_mays.AGPv4.dna.toplevel_sorted.chrom.sizes \
# > ${outFolder}/${sample_prefix}_Olap_gene_space_V.bed
# 
# cd ${outFolder}
# 
# bedtools closest \
# -a $DE_UMR_bed  \
# -b ~/refseqs/maize/B73.structuralTEv2.2018.12.20.filteredTE.disjoined_sup_sorted_bed6.bed \
# -mdb all \
# -t all \
# -D b \
# -g ~/refseqs/maize/Zea_mays.AGPv4.dna.toplevel_sorted.chrom.sizes \
# > ${outFolder}/${sample_prefix}_Olap_TE_order.bed
# 
# module load R/3.5.0-gnu
# 
# # load R
# R
# 
# # data
# path_to_data_files="B_M_bams_2"
# contrast = "B73_leaf3_V3.vs.Mo17_leaf3_V3"
# filter_FC = 3
# bin_size = 100
# window_spacing = 50
###########################
########
###########################

library(tidyverse)
library(ggthemes)
# library("seqinr")
old.scipen <- getOption("scipen")
options(scipen=999)
library(wesanderson)
library(RColorBrewer)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8))

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

###### set up
projectFolder <- paste0(path_to_data_files, "_CSAW")
dir.create(projectFolder)

data_dir <- paste0(projectFolder, "/", contrast, "_FCF", filter_FC, "_bin", bin_size, "_space", window_spacing)

outFolder <- paste0(projectFolder, "/", contrast, "_FCF", filter_FC, "_bin", bin_size, "_space", window_spacing, "/annotation")
dir.create(outFolder)


# Gene ONLY proximity (TE afterwards) annotation

# Only gene are overlapped in this version to get an accurate gene-distal list

# Make three categories

# 1. genic
# 2. gene_proximal < 2kb (ths is captured by 1 or 2kb up/down categories)
# 3. gene distal > 2kb (anything called intergenic)

# Header:
# Column1: chromosome
# Column2: start
# Column3: stop
# Column4: UMR_genotype
# Column5: log2FC - log2 fold change of best window in the DE-UMR region
# Column6: strand

######################
## overlaps
######################

#### resolve multi-overlaps and summarise

sample_to_crunch = contrast
reference_overlap = "gene_space_V"
#unmeth_Olap_Hv-rep1-H2-ACRs.bed

overlaps <- read_tsv(paste0(outFolder, "/", sample_to_crunch, "_Olap_", reference_overlap, ".bed"),
               col_names = c("chr", "start", "end", "UMR_genotype", "log2FC", "strand", 
                               "b_chr", "b_start", "b_end", "b_name", "b_score", "b_strand", "distance"), 
               cols(chr = col_character(),
                    # Olap_file = col_character(), 
                    b_chr = col_character()))

overlaps 
# 51,293

# resolve multi-overlaps
# I think most of the multi-overlaps will be at distance zero - this pipeline is mainly aimed at resolving that. Its possible some will be due to a peak being exactly the same length between an upstream and downstream gene, but I might just pick one for these because I think they will be rare
# get distinct tile + feature rows
overlaps_distinct <- overlaps %>%
  select(chr, start, end, b_name, distance, UMR_genotype, log2FC) %>%
  distinct()
# 50,934 

overlaps %>%
  select(chr, start, end) %>%
  distinct()
# 50,618


# collapse, I take mean distance for tiles with multioverlaps. I think this is ok because gene-gene and TE-TE will called ambiguous because they are half way and for gene-TE if its distance is 0 then I will call it geneic TE and mean(0,0 = 0) and if distance is > 0 I'll just call it ambiguous. This scheme ensures no tile is double counted for distribution purposes.
overlaps_distinct_collapsed <- overlaps_distinct %>%
  group_by(chr, start, end, UMR_genotype, log2FC) %>%
  arrange(chr, start, b_name) %>%
  summarise(feature = paste(b_name, collapse = "-"), distance2 = mean(distance))

overlaps_distinct_collapsed
# 50,618

overlaps_distinct_collapsed %>% group_by(feature) %>% summarise(n = n()) %>% print(n = 35)

# # A tibble: 19 × 2
#    feature                                         n
#    <chr>                                       <int>
#  1 .                                             121
#  2 lincRNA_gene                                 2469
#  3 lincRNA_gene-lincRNA_gene                       1
#  4 lincRNA_gene-miRNA_gene                         1
#  5 lincRNA_gene-nonSyntenic_gene                  48
#  6 lincRNA_gene-nonSyntenic_gene-syntenic_gene     1
#  7 lincRNA_gene-syntenic_gene                     63
#  8 lincRNA_gene-tRNA_gene                          1
#  9 miRNA_gene                                     71
# 10 miRNA_gene-miRNA_gene                          70
# 11 miRNA_gene-nonSyntenic_gene                     1
# 12 miRNA_gene-syntenic_gene                        1
# 13 nonSyntenic_gene                            16409
# 14 nonSyntenic_gene-syntenic_gene                101
# 15 nonSyntenic_gene-syntenic_gene-tRNA_gene        1
# 16 nonSyntenic_gene-tRNA_gene                      4
# 17 syntenic_gene                               30180
# 18 syntenic_gene-tRNA_gene                        21
# 19 tRNA_gene                                    1054

############ Gene overlaps 


# sRNA data annotation rules
# miRNA > syntenic_gene > non_syntenic > TE > ncRNA

# This distance does now take into account gene strand...

overlaps_distinct_collapsed_filtered <- overlaps_distinct_collapsed %>% ungroup() %>% 
  # slice(1:10000) %>%
  mutate(classification = ifelse(grepl("miRNA_gene", feature) & distance2 == 0, "miRNA_gene",
                          ifelse(grepl("syntenic_gene", feature) & distance2 == 0, "syntenic_gene",
                          ifelse(grepl("nonSyntenic_gene", feature) & distance2 == 0, "nonSyntenic_gene",
                          ifelse(grepl("tRNA_gene", feature) & distance2 == 0, "tRNA_gene",
                          ifelse(grepl("lincRNA_gene", feature) & distance2 == 0, "lincRNA_gene",
                          ifelse(grepl("_gene", feature) & distance2 < 0 & distance2 >= -1000, "1kb_upstream_gene",
                          ifelse(grepl("_gene", feature) & distance2 < -1000 & distance2 >= -2000, "2kb_upstream_gene",
                          ifelse(grepl("_gene", feature) & distance2 > 0 & distance2 <= 1000, "1kb_downstream_gene",
                          ifelse(grepl("_gene", feature) & distance2 > 1000 & distance2 <= 2000, "2kb_downstream_gene", "intergenic"))))))))))

overlaps_distinct_collapsed_filtered
# 50,618

#######
# write 100 bp annotation?
overlaps_distinct_collapsed_filtered_out <- overlaps_distinct_collapsed_filtered %>%
  select(chr, start, end, UMR_genotype, log2FC, classification)

write.table(overlaps_distinct_collapsed_filtered_out, paste0(outFolder,"/",sample_to_crunch, "_", reference_overlap, "_annotated.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

########
#summarise
overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed_filtered %>%
  group_by(classification) %>%
  summarise(n = n())

overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed_filtered_summary %>% 
  mutate(percentage = n/sum(n)*100,
         Mb = n * 100 / 1000000)

overlaps_distinct_collapsed_filtered_summary
# # A tibble: 10 × 4
#    classification          n percentage     Mb
#    <chr>               <int>      <dbl>  <dbl>
#  1 1kb_downstream_gene  2617      5.17  0.262 
#  2 1kb_upstream_gene    2370      4.68  0.237 
#  3 2kb_downstream_gene  1579      3.12  0.158 
#  4 2kb_upstream_gene    1781      3.52  0.178 
#  5 intergenic          27996     55.3   2.80  
#  6 lincRNA_gene          478      0.944 0.0478
#  7 miRNA_gene             81      0.160 0.0081
#  8 nonSyntenic_gene     3867      7.64  0.387 
#  9 syntenic_gene        9775     19.3   0.978 
# 10 tRNA_gene              74      0.146 0.0074


write_tsv(overlaps_distinct_collapsed_filtered_summary, paste0(outFolder, "/", sample_to_crunch, "_olap_gene_summary_2.tsv"))

####### pull distal

distal_only_bed <- overlaps_distinct_collapsed_filtered_out %>% filter(classification == "intergenic") %>% 
  filter(chr %in% c(1:10)) %>% # no contigs
  mutate(name = sample_to_crunch) %>% select(chr, start, end, UMR_genotype, log2FC, classification)
# 27,687

# write
write.table(distal_only_bed, paste0(outFolder, "/", sample_to_crunch, "_distal.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

# If i want to look at TEs now I could merge in the TE order annotation as above to work out the percent of distal UMRs that fall in TEs

### Add TE family annotation
# In this version I merge in TE order/superfamily annotations from the disjoined file.
# This version include miRNAs in B73 refgen and miRNAs in miRbase and other annotated non-coding RNAs like linc RNAs and t-RNAs 
# This version also calssifies genes as syntenic or non syntenic

#### overlaps

# Just overlap with the TE reference for define the TE for each tile then merge this into IV above to annotate TE type.

#### summarise

# resolve multi-overlaps
reference_overlap = "TE_order"

overlaps <- read_tsv(paste0(outFolder, "/", sample_to_crunch, "_Olap_", reference_overlap, ".bed"),
               col_names = c("chr", "start", "end", "UMR_genotype", "log2FC", "strand",
                               "b_chr", "b_start", "b_end", "b_name", "b_score", "b_strand", "distance"), 
               cols(chr = col_character(),
                    b_chr = col_character()) 
               # , n_max = 100000
               )

overlaps 
# 52,275

overlaps <- overlaps %>% mutate(order = substr(b_name, start = 1, stop = 2))

print(overlaps, width = Inf)

# resovle multi-overlaps
# I think most of the multi-overlaps will be at distance zero - this pipeline is mainly aimed at resolving that. Its possible some will be due to a peak being exactally the same length between an upstrea and downstream gene, but I might just pick one for these because I think they will be rare
# get distinct tile + feature rows
overlaps_distinct <- overlaps %>%
  select(chr, start, end, order, distance, UMR_genotype, log2FC) %>%
  distinct()
# 51,646

# colapse, I take mean distance for tiles with multioverlaps. I think this is ok because gene-gene and TE-TE will called ambiguous becaus ethey are half way and for gene-TE if its distance is 0 then I will call it geneic TE and mean(0,0 = 0) and if distance is > 0 I'll just call it ambiguous. This scheme ensures no tile is double counted for distribution purposes.
overlaps_distinct_collapsed <- overlaps_distinct %>%
  group_by(chr, start, end, UMR_genotype, log2FC) %>%
  # arrange(chr, start, b_name) %>%
  summarise(feature = paste(order, collapse = "-"), distance2 = mean(distance))

overlaps_distinct_collapsed
# 350,618

overlaps_distinct_collapsed %>% group_by(feature) %>% summarise(n = n()) %>% print(n = 123)
# # A tibble: 21 × 2
#    feature      n
#    <chr>    <int>
#  1 DH        5229
#  2 DH-DT       56
#  3 DH-RL       49
#  4 DH-RL-DT     1
#  5 DH-RS        3
#  6 DT       16397
#  7 DT-DH       50
#  8 DT-DH-RL     2
#  9 DT-DT        1
# 10 DT-RL      305
# 11 DT-RS        1
# 12 RI         279
# 13 RI-DT        1
# 14 RL       27521
# 15 RL-DH       54
# 16 RL-DT      492
# 17 RL-RI        1
# 18 RL-RS        4
# 19 RS         167
# 20 RS-DT        3
# 21 RS-RL        2

############ Gene and TE overlaps 


# sRNA data annotation rules
# miRNA > syntenic_gene > non_syntenic > TE > ncRNA

overlaps_distinct_collapsed_filtered <- overlaps_distinct_collapsed %>% ungroup() %>% 
  # slice(1:10000) %>%
  mutate(classification = ifelse(nchar(feature) == 2 & distance2 == 0, feature, 
                          ifelse(distance2 == 0, "multi_TE", NA)))

overlaps_distinct_collapsed_filtered

#######
# write 100 bp annotation?
overlaps_distinct_collapsed_filtered_out <- overlaps_distinct_collapsed_filtered %>%
  select(chr, start, end, UMR_genotype, log2FC, classification)

write.table(overlaps_distinct_collapsed_filtered_out, paste0(outFolder,"/",sample_to_crunch, "_", reference_overlap,  "_annotated.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

# write a version with NA tiles removed
overlaps_distinct_collapsed_filtered_out_no_na <- overlaps_distinct_collapsed_filtered_out %>% na.omit()
overlaps_distinct_collapsed_filtered_out_no_na
# # A tibble: 19,545
# 70% smaller
write.table(overlaps_distinct_collapsed_filtered_out_no_na, paste0(outFolder,"/",sample_to_crunch,"_", reference_overlap, "_annotated_no_na.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

########
#summarise
overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed_filtered %>%
  group_by(classification) %>%
  summarise(n = n())

overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed_filtered_summary %>% 
  mutate(percentage = n/sum(n)*100,
         Mb = n * 100 / 1000000)

overlaps_distinct_collapsed_filtered_summary

# # A tibble: 7 × 4
#   classification     n percentage     Mb
#   <chr>          <int>      <dbl>  <dbl>
# 1 DH              2599     5.13   0.260 
# 2 DT              2203     4.35   0.220 
# 3 multi_TE        1025     2.02   0.102 
# 4 RI                51     0.101  0.0051
# 5 RL             13641    26.9    1.36  
# 6 RS                26     0.0514 0.0026
# 7 NA             31073    61.4    3.11   

write_tsv(overlaps_distinct_collapsed_filtered_summary, paste0(outFolder, "/", sample_to_crunch, "_", reference_overlap,"_summary_2.tsv"))

#######################################
#### MERGE INTO EXISTING gene OVERLAP
#######################################

# TE ref
reference_overlap2 = "TE_order"

TE_ref <- read_tsv(paste0(outFolder,"/",sample_to_crunch, "_", reference_overlap2, "_annotated_no_na.bed"), col_names = c("chr", "start", "end", "UMR_genotype", "log2FC", "classification"),
                   cols(chr = col_character()))

# gene ref
reference_overlap = "gene_space_V"

ref <- read_tsv(paste0(outFolder,"/",sample_to_crunch, "_", reference_overlap, "_annotated.bed"), col_names = c("chr", "start", "end", "UMR_genotype", "log2FC", "feature"),
                cols(chr = col_character()))

# TE_ref <- overlaps_distinct_collapsed_filtered_out_no_na
# or


# TE_ref <- TE_ref %>% rename(classification = name)

# TE_ref <- TE_ref %>% select(-score, -strand)

ref_TE_order <- ref %>% left_join(TE_ref, by = c("chr", "start", "end", "UMR_genotype", "log2FC"))
# 50,618

ref_TE_super_compare <- ref_TE_order %>% mutate(joins = paste0(feature, "_", classification))

ref_TE_super_compare_summary <- ref_TE_super_compare %>% 
  group_by(joins) %>%
  summarise(n = n())

# hmm this suggests that there are intergenic loci that actually do have TEs (this shouldnt be the case)? - I cant remember why I thought this anymore
print(ref_TE_super_compare_summary, n = 59)

ref_TE_order_out <- ref_TE_order %>% dplyr::rename(TE = classification)

write.table(ref_TE_order_out, paste0(outFolder,"/", sample_to_crunch, "_", reference_overlap, "_", reference_overlap2, ".bed"), sep = "\t", quote = F, row.names = F, col.names = F)

ref_TE_order_out %>% group_by(feature) %>% summarise(n = n())

#### Add meta annotation (local)

# NOTE: THIS CURRENTLY CALLS AN ACR GENIC IF IT OVERLAPS THE GENE. FOR UMRs WE CALLED ANY OVERLAP OF THE PROMOTER A PROXIMAL NOT GENIC... TO PROPERLY COMPARE I SHOUDL RECALL BOTH USING THE SAME SCHEME

features <- ref_TE_order_out

features_meta <- features %>% 
  mutate(meta_feature = ifelse(feature %in% c("1kb_upstream_gene", "2kb_upstream_gene", "1kb_downstream_gene", "2kb_downstream_gene"), "gene_proximal", 
                               ifelse(feature %in% c("syntenic_gene", "nonSyntenic_gene"), "protein_coding_gene", 
                                      ifelse(feature %in% c("miRNA_gene", "tRNA_gene", "lincRNA_gene"), "non_coding_RNA", 
                                             ifelse(feature %in% c("intergenic"), "gene_distal", feature)))))

features_meta %>% distinct(meta_feature)

# # A tibble: 4 × 1
#   meta_feature       
#   <chr>              
# 1 gene_distal        
# 2 protein_coding_gene
# 3 gene_proximal      
# 4 non_coding_RNA    

features_meta %>% distinct(feature)

features_meta %>% distinct(feature, meta_feature)

# # A tibble: 10 × 2
#    feature             meta_feature       
#    <chr>               <chr>              
#  1 intergenic          gene_distal        
#  2 nonSyntenic_gene    protein_coding_gene
#  3 2kb_upstream_gene   gene_proximal      
#  4 1kb_upstream_gene   gene_proximal      
#  5 1kb_downstream_gene gene_proximal      
#  6 syntenic_gene       protein_coding_gene
#  7 2kb_downstream_gene gene_proximal      
#  8 lincRNA_gene        non_coding_RNA     
#  9 miRNA_gene          non_coding_RNA     
# 10 tRNA_gene           non_coding_RNA         

write.table(features_meta, paste0(outFolder,"/", sample_to_crunch, "_", reference_overlap, "_", reference_overlap2, "_meta_feature.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

######################
#############
# subset - distal only
features_meta %>% group_by(meta_feature) %>% summarise(n = n()) %>% print(n = 30)

distal_simple <- features_meta %>%
   mutate(classification = ifelse(grepl(paste("gene_proximal", "non_coding_RNA", "protein_coding_gene", sep = "|"), meta_feature), "gene_proximal", "gene_distal"))
# mutate(classification = ifelse(grepl("protein_coding_gene", Mfeatures), "gene_proximal", "gene_distal"))

meta_feature_summay <- distal_simple %>% 
  filter(chr %in% c(1:10)) %>% # no contigs
  group_by(meta_feature) %>%
  summarise(n = n()) %>%
  mutate(percent = n/sum(n)*100)

meta_feature_summay

# # A tibble: 4 × 3
#   meta_feature            n percent
#   <chr>               <int>   <dbl>
# 1 gene_distal         27687   55.2 
# 2 gene_proximal        8307   16.6 
# 3 non_coding_RNA        633    1.26
# 4 protein_coding_gene 13542   27.0 

write.table(meta_feature_summay, paste0(outFolder, "/", sample_to_crunch, "_distal_summary.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

############
# summarise features
############

# genes
features_meta %>% 
  mutate(domain = "DE_UMRs") %>%
  group_by(domain, meta_feature) %>%
  summarise(tiles = n()) %>%
  group_by(domain) %>%
  mutate(total_tiles = sum(tiles),
         percent = tiles/total_tiles*100)

# # A tibble: 4 × 5
# # Groups:   domain [1]
#   domain  meta_feature        tiles total_tiles percent
#   <chr>   <chr>               <int>       <int>   <dbl>
# 1 DE_UMRs gene_distal         27996       50618   55.3 
# 2 DE_UMRs gene_proximal        8347       50618   16.5 
# 3 DE_UMRs non_coding_RNA        633       50618    1.25
# 4 DE_UMRs protein_coding_gene 13642       50618   27.0 

# genes
features_meta %>% 
  mutate(domain = "DE_UMRs") %>%
  group_by(UMR_genotype, meta_feature) %>%
  summarise(tiles = n()) %>%
  group_by(UMR_genotype) %>%
  mutate(total_tiles = sum(tiles),
         percent = tiles/total_tiles*100)

# # A tibble: 12 × 5
# # Groups:   UMR_genotype [3]
#    UMR_genotype  meta_feature        tiles total_tiles percent
#    <chr>         <chr>               <int>       <int>   <dbl>
#  1 B73_leaf3_V3  gene_distal         20433       36008  56.7  
#  2 B73_leaf3_V3  gene_proximal        6016       36008  16.7  
#  3 B73_leaf3_V3  non_coding_RNA        550       36008   1.53 
#  4 B73_leaf3_V3  protein_coding_gene  9009       36008  25.0  
#  5 mixed         gene_distal            17          85  20    
#  6 mixed         gene_proximal          10          85  11.8  
#  7 mixed         non_coding_RNA          2          85   2.35 
#  8 mixed         protein_coding_gene    56          85  65.9  
#  9 Mo17_leaf3_V3 gene_distal          7546       14525  52.0  
# 10 Mo17_leaf3_V3 gene_proximal        2321       14525  16.0  
# 11 Mo17_leaf3_V3 non_coding_RNA         81       14525   0.558
# 12 Mo17_leaf3_V3 protein_coding_gene  4577       14525  31.5  

# genes
mC_domains2_anno_summary <- features_meta %>% 
  mutate(domain = "DE_UMRs") %>%
  group_by(domain, feature) %>%
  summarise(tiles = n(), meta_feature = unique(meta_feature)) %>%
  group_by(domain) %>%
  mutate(total_tiles = sum(tiles),
         percent = tiles/total_tiles*100)
  
print(mC_domains2_anno_summary, n = 40)
# # A tibble: 10 × 6
# # Groups:   domain [1]
#    domain  feature             tiles meta_feature        total_tiles percent
#    <chr>   <chr>               <int> <chr>                     <int>   <dbl>
#  1 DE_UMRs 1kb_downstream_gene  2617 gene_proximal             50618   5.17 
#  2 DE_UMRs 1kb_upstream_gene    2370 gene_proximal             50618   4.68 
#  3 DE_UMRs 2kb_downstream_gene  1579 gene_proximal             50618   3.12 
#  4 DE_UMRs 2kb_upstream_gene    1781 gene_proximal             50618   3.52 
#  5 DE_UMRs intergenic          27996 gene_distal               50618  55.3  
#  6 DE_UMRs lincRNA_gene          478 non_coding_RNA            50618   0.944
#  7 DE_UMRs miRNA_gene             81 non_coding_RNA            50618   0.160
#  8 DE_UMRs nonSyntenic_gene     3867 protein_coding_gene       50618   7.64 
#  9 DE_UMRs syntenic_gene        9775 protein_coding_gene       50618  19.3  
# 10 DE_UMRs tRNA_gene              74 non_coding_RNA            50618   0.146

mC_domains2_anno_summary %>% ungroup() %>% filter(is.na(feature)) %>% summarise(total = sum(tiles))
# there are zero!

write.table(mC_domains2_anno_summary, paste0(outFolder,"/", sample_to_crunch, "_domains_features", 
                                    "_freq.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

# TEs

mC_domains2_anno_summary_TE <-  features_meta %>% 
  mutate(domain = "DE_UMRs") %>%
  mutate(TE = ifelse(is.na(TE), "non-TE", TE)) %>%
  group_by(domain, TE) %>%
  summarise(tiles = n()) %>%
  group_by(domain) %>%
  mutate(total_tiles = sum(tiles),
         percent = tiles/total_tiles*100)
  
write.table(mC_domains2_anno_summary_TE, paste0(outFolder,"/", sample_to_crunch, "_TE_features", 
                                    "_freq.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

# # A tibble: 7 × 5
# # Groups:   domain [1]
#   domain  TE       tiles total_tiles percent
#   <chr>   <chr>    <int>       <int>   <dbl>
# 1 DE_UMRs DH        2599       50618  5.13  
# 2 DE_UMRs DT        2203       50618  4.35  
# 3 DE_UMRs multi_TE  1025       50618  2.02  
# 4 DE_UMRs non-TE   31073       50618 61.4   
# 5 DE_UMRs RI          51       50618  0.101 
# 6 DE_UMRs RL       13641       50618 26.9   
# 7 DE_UMRs RS          26       50618  0.0514

# distal TEs
mC_domains2_anno_summary_TE_distal <-  features_meta %>% 
  mutate(domain = "DE_UMRs") %>%
  mutate(TE = ifelse(is.na(TE), "non-TE", TE)) %>%
  filter(feature == "intergenic") %>%
  group_by(domain, TE) %>%
  summarise(tiles = n()) %>%
  group_by(domain) %>%
  mutate(total_tiles = sum(tiles),
         percent = tiles/total_tiles*100)
  
print(mC_domains2_anno_summary_TE, n = 40)

write.table(mC_domains2_anno_summary_TE_distal, paste0(outFolder,"/", sample_to_crunch, "_distal_TE_features", 
                                    "_freq.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)


#####plots gene features

plot_colours_bar <- c( "#B2182B", "#EF8A62", "lightgrey", "black", "#67A9CF", "#2166AC", "green" , "purple", "red", "darkgrey")
plot_order_bar <- c("2kb_upstream_gene", "1kb_upstream_gene", "nonSyntenic_gene", "syntenic_gene", "1kb_downstream_gene", "2kb_downstream_gene", "lincRNA_gene",   "tRNA_gene", "miRNA_gene", "intergenic"
                    #, "ncRNA_gene", "miRNA_gene",  "TE_DH", "TE_DT", "TE_RI", "TE_RL", "TE_RS", "TE_other"
                    )
# domain_order <- c("Heterochromatin", "RdDM", "CG_only", "Unmethylated", "Intermediate", "no_sites", "Missing_Data")

mC_domains2_anno_summary_order <- mC_domains2_anno_summary %>% ungroup() %>%
  filter(!is.na(feature)) %>% # remove the 2K tiles with no features
  # mutate(feature = ifelse(feature %in% c("TE_NA", "TE_multi_TE"), "TE_other", feature))  %>% # remove the TE tiles with no features
  mutate(feature = factor(feature, levels = plot_order_bar))

# stacked bar (no eror bars) scales free
g1 <- ggplot(mC_domains2_anno_summary_order, aes(x = domain, y = percent, fill = feature)) +
  geom_bar(stat = "identity") +
  # geom_errorbar(aes(ymin=average_percent-sd_normalised, ymax=average_percent+sd_normalised), width=.3) +
  # facet_grid(tissue ~ .) +
  scale_fill_manual(values = plot_colours_bar, breaks = plot_order_bar) +
  # facet_grid (meta_feature ~ ., scales = 'free') +
  theme_minimal() +
  text_size_theme_8
# print(g1)

# g1

ggsave(filename = paste0(outFolder, "/", sample_to_crunch, "_feature_freq_domain_stack_bar_no_TE_scales_fixed.pdf"), plot = g1, h = 3, w = 2.7)

#####plots TEs features

plot_colours_bar <- c("#67A9CF", wes_palette("IsleofDogs1", 6, type = "discrete"))

plot_order_bar <- c("non-TE", "DH", "DT", "RI", "RL", "RS", "multi_TE")

mC_domains2_anno_summary_TE_order <- mC_domains2_anno_summary_TE %>% ungroup() %>%
  filter(!is.na(TE)) %>% # remove the 2K tiles with no features
  # mutate(feature = ifelse(feature %in% c("TE_NA", "TE_multi_TE"), "TE_other", feature))  %>% # remove the TE tiles with no features
  mutate(TE = factor(TE, levels = plot_order_bar))


# stacked bar (no eror bars) scales free
g1 <- ggplot(mC_domains2_anno_summary_TE_order, aes(x = domain, y = percent, fill = TE)) +
  geom_bar(stat = "identity") +
  # geom_errorbar(aes(ymin=average_percent-sd_normalised, ymax=average_percent+sd_normalised), width=.3) +
  # facet_grid(tissue ~ .) +
  scale_fill_manual(values = plot_colours_bar, breaks = plot_order_bar) +
  # facet_grid (meta_feature ~ ., scales = 'free') +
  theme_minimal() +
  text_size_theme_8
# print(g1)

# g1

ggsave(filename = paste0(outFolder, "/", sample_to_crunch, "_feature_freq_domain_stack_bar_TE_scales_fixed.pdf"), plot = g1, h = 3, w = 1.8)

#####plots distal TE features

plot_colours_bar <- c("#67A9CF", wes_palette("IsleofDogs1", 6, type = "discrete"))

plot_order_bar <- c("non-TE", "DH", "DT", "RI", "RL", "RS", "multi_TE")

mC_domains2_anno_summary_TE_order <- mC_domains2_anno_summary_TE_distal %>% ungroup() %>%
  filter(!is.na(TE)) %>% # remove the 2K tiles with no features
  # mutate(feature = ifelse(feature %in% c("TE_NA", "TE_multi_TE"), "TE_other", feature))  %>% # remove the TE tiles with no features
  mutate(TE = factor(TE, levels = plot_order_bar))


# stacked bar (no eror bars) scales free
g1 <- ggplot(mC_domains2_anno_summary_TE_order, aes(x = domain, y = percent, fill = TE)) +
  geom_bar(stat = "identity") +
  # geom_errorbar(aes(ymin=average_percent-sd_normalised, ymax=average_percent+sd_normalised), width=.3) +
  # facet_grid(tissue ~ .) +
  scale_fill_manual(values = plot_colours_bar, breaks = plot_order_bar) +
  # facet_grid (meta_feature ~ ., scales = 'free') +
  theme_minimal() +
  text_size_theme_8
# print(g1)

# g1

ggsave(filename = paste0(outFolder, "/", sample_to_crunch, "_feature_freq_domain_stack_bar_distal_TE_scales_fixed.pdf"), plot = g1, h = 3, w = 1.8)


##### plots all and distal TE features

combo_plot <- mC_domains2_anno_summary_TE_distal %>% ungroup() %>% mutate(domain = "distal_DE_UMRs")  %>% bind_rows(ungroup(mC_domains2_anno_summary_TE))


plot_colours_bar <- c("#67A9CF", wes_palette("IsleofDogs1", 6, type = "discrete"))

plot_order_bar <- c("non-TE", "DH", "DT", "RI", "RL", "RS", "multi_TE")

mC_domains2_anno_summary_TE_order <- combo_plot %>% ungroup() %>%
  filter(!is.na(TE)) %>% # remove the 2K tiles with no features
  # mutate(feature = ifelse(feature %in% c("TE_NA", "TE_multi_TE"), "TE_other", feature))  %>% # remove the TE tiles with no features
  mutate(TE = factor(TE, levels = plot_order_bar))


# stacked bar (no eror bars) scales free
g1 <- ggplot(mC_domains2_anno_summary_TE_order, aes(x = domain, y = percent, fill = TE)) +
  geom_bar(stat = "identity") +
  # geom_errorbar(aes(ymin=average_percent-sd_normalised, ymax=average_percent+sd_normalised), width=.3) +
  # facet_grid(tissue ~ .) +
  scale_fill_manual(values = plot_colours_bar, breaks = plot_order_bar) +
  # facet_grid (meta_feature ~ ., scales = 'free') +
  theme_minimal() +
  text_size_theme_8
# print(g1)

# g1

ggsave(filename = paste0(outFolder, "/", sample_to_crunch, "_feature_freq_domain_stack_bar_distal_combo_TE_scales_fixed.pdf"), plot = g1, h = 3, w = 2.0)

