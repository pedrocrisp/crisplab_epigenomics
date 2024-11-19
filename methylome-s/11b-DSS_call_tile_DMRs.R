#!/usr/bin/Rscript
##########

# Peter Crisp
# 2018-07-09
# R script to call DMRs using DSS, supplying 100 bp tiles

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

####### libs
library(tidyverse)
library(DSS)
library(bsseq)
library(limma)
library(tidygenomics)
library(ggthemes)
old.scipen <- getOption("scipen")
options(scipen=999)

#######  args

# contrast = "PAC004_root_HC_mother.vs.PAC012_gen_ST_0"
# DMR_contrasts_table_file = "DMR_tests_combos_all_table.tsv"
# path_to_data_files = "analysis/tiles_filtered_4C_2x"

####### set up

outFolder <- paste0(path_to_data_files, "_DSS_DMRs")
dir.create(outFolder)

DMR_contrasts_table <- read_tsv(DMR_contrasts_table_file)

sample1 = pull(filter(DMR_contrasts_table, test_name == contrast), sample1)
sample1

sample2 = pull(filter(DMR_contrasts_table, test_name == contrast), sample2)
sample2

print(paste0("Calling DMR tiles for ", sample1, " vs ", sample2))

#######  #######
# CG
#######  #######

context = "CG"

#######  read in the files
sample1_data <- read_tsv(file.path(path_to_data_files, paste0(sample1, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample1_data <- sample1_data %>% select("chr", "pos", "N", "X", "sites")
# sample1_data <- sample1_data %>% filter(chr == "Chr01")
sample1_data

# # A tibble: 3,456,028 × 5
#    chr     pos     N     X sites
#    <chr> <dbl> <dbl> <dbl> <dbl>
#  1 Chr01   401    34    28     5
#  2 Chr01   601    61    36     6
#  3 Chr01   901    89     0     8
#  4 Chr01  1001    97     1     6
#  5 Chr01  1101    73    59     4

sample2_data <- read_tsv(file.path(path_to_data_files, paste0(sample2, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample2_data <- sample2_data %>% select("chr", "pos", "N", "X", "sites")
# sample2_data <- sample2_data %>% filter(chr == "Chr01")
sample2_data

# # A tibble: 3,340,994 × 5
#    chr     pos     N     X sites
#    <chr> <dbl> <dbl> <dbl> <dbl>
#  1 Chr01   601    57    40     6
#  2 Chr01   901    70     0     8
#  3 Chr01  1001    71     2     6
#  4 Chr01  1101    46    33     4
#  5 Chr01  1301    51    41     4

#######  1. create an object of BSseq class (requires bsseq Bioconductor package)

BSobj <- makeBSseqData(list(sample1_data, sample2_data), c(sample1, sample2) )
BSobj

# An object of type 'BSseq' with
#   3500080 methylation loci
#   2 samples
# has not been smoothed
# All assays are in-memory

####### 2 Perform statistical test for DML by calling DMLtest function
# cant remember why I didnt use smoothing? I guess because they are 100 bp tiles,
# we have already smoothed into these tiles
# t1 <- proc.time()
dmlTest <- DMLtest(BSobj, group1=sample1, group2=sample2, smoothing = F, equal.disp = T, ncores = 4)
# proc.time() - t1

# user  system elapsed
# 508.827   0.139 509.181

####### 3. Call DMRs

dmls <- callDML(dmlTest, p.threshold=0.01, delta = 0.1)
dmls_tbl <- as.tibble(dmls)
dmls_tbl

# # A tibble: 16,163 × 12
#    chr        pos     mu1    mu2   diff diff.se   stat    phi1    phi2     pval
#    <chr>    <int>   <dbl>  <dbl>  <dbl>   <dbl>  <dbl>   <dbl>   <dbl>    <dbl>
#  1 Chr01  4406901 0.762   0.0844  0.678  0.0644  10.5  0.00674 0.00674 7.09e-26
#  2 Chr01 20805301 0.882   0.215   0.667  0.0641  10.4  0.00674 0.00674 2.34e-25
#  3 Chr01 22843701 0.0245  0.793  -0.768  0.0629 -12.2  0.00674 0.00674 2.65e-34
#  4 Chr01 24817401 0.914   0.170   0.745  0.0545  13.7  0.00674 0.00674 1.60e-42
#  5 Chr01 57308101 0.0447  0.708  -0.663  0.0640 -10.4  0.00674 0.00674 3.58e-25

# merge adjacent DMRs
# note: the bed must be ordered - tidy genomics does not check!
dmls_bed <- dmls_tbl %>% mutate(chr = as.character(chr), start = as.double(pos), end = pos+100, type = ifelse(diff > 0, "hyper", "hypo")) %>%
  select(chr, start, end, diff, type)  %>%
  arrange(chr, start)

dmls_bed
# # A tibble: 16,163 × 5
#    chr    start    end   diff type 
#    <chr>  <dbl>  <dbl>  <dbl> <chr>
#  1 Chr01   5701   5801  0.408 hyper
#  2 Chr01  35901  36001 -0.440 hypo 
#  3 Chr01  62101  62201  0.344 hyper
#  4 Chr01  78401  78501 -0.354 hypo 
#  5 Chr01  87301  87401  0.330 hyper

dmls_bed_summary <- dmls_bed %>%
  group_by(type) %>%
  mutate(length = end - start) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(diff))
dmls_bed_summary

# # A tibble: 2 × 5
#   type      n mean_length max_length mean_diff
#   <chr> <int>       <dbl>      <dbl>     <dbl>
# 1 hyper  8613         100        100     0.362
# 2 hypo   7550         100        100    -0.357

dmls_bed_hyper <- filter(dmls_bed, type == "hyper")
dmls_bed_hypo <- filter(dmls_bed, type == "hypo")

dmls_bed_hyper <- genome_cluster(dmls_bed_hyper, by=c("chr", "start", "end"), max_distance = 1)
dmls_bed_hypo <- genome_cluster(dmls_bed_hypo, by=c("chr", "start", "end"), max_distance = 1)

dmls_bed_pre_merge <- bind_rows(dmls_bed_hyper, dmls_bed_hypo)
# A tibble: 16,163 × 6

# dmls_bed %>% select(chr, start, end, diff) %>% arrange(chr, start)
# dmls_bed %>% select(chr, start, end, diff, cluster_id) %>% arrange(cluster_id)

dmls_bed_merged <- dmls_bed_pre_merge %>% group_by(cluster_id, type) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), mean_diff = mean(diff), max_diff = max(diff), min_diff = min(diff)) %>%
  mutate(length = end - start) %>%
  ungroup()
dmls_bed_merged

# # A tibble: 15,479 × 9
#    cluster_id type  chr      start      end mean_diff max_diff min_diff length
#         <dbl> <chr> <chr>    <dbl>    <dbl>     <dbl>    <dbl>    <dbl>  <dbl>
#  1          0 hyper Chr01     5701     5801     0.408    0.408    0.408    100
#  2          0 hypo  Chr01    35901    36001    -0.440   -0.440   -0.440    100
#  3          1 hyper Chr01    62101    62201     0.344    0.344    0.344    100
#  4          1 hypo  Chr01    78401    78501    -0.354   -0.354   -0.354    100
#  5          2 hyper Chr01   725901   726001     0.350    0.350    0.350    100

# dmls_bed_merged %>% filter(length > 100)

# no directionality to DMR 59,731; if hypo/hyper is considered 60,134

dmls_bed_merge_summarised <- dmls_bed_merged %>%
  group_by(type) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(mean_diff))
dmls_bed_merge_summarised

# # A tibble: 2 × 5
#   type      n mean_length max_length mean_diff
#   <chr> <int>       <dbl>      <dbl>     <dbl>
# 1 hyper  8244        104.        800     0.361
# 2 hypo   7235        104.        500    -0.356

####### 4. Make results files and write out tables

dmlTest_calls <- as.tibble(dmlTest) %>% left_join(dmls_tbl, by = colnames(dmlTest))
dmlTest_calls

### write DMR calls (large) file - only retain mu1 and mu2 cut off their columns to make file smaller
dmlTest_calls_out <- dmlTest_calls %>% select(chr, pos, mu1, mu2)
write_tsv(dmlTest_calls_out, paste0(outFolder, "/", contrast, "_", context, "_dmlTest_calls.tsv"))

### write only DMR tiles
write_csv(dmls_tbl, paste0(outFolder, "/", contrast, "_", context, "_DMRs.csv"))

### write summary of DMR tiles
write_csv(dmls_bed_summary, paste0(outFolder, "/", contrast, "_", context, "_DMR_tiles_summary.csv"))

### write merged DMRs
dmls_bed_merged_out <- dmls_bed_merged %>% select(chr, start, end, type, mean_diff, max_diff, min_diff, length)
write_tsv(dmls_bed_merged_out, paste0(outFolder, "/", contrast, "_", context, "_DMRs_merged.tsv"))

### write summary of merged DMR tiles
write_csv(dmls_bed_merge_summarised, paste0(outFolder, "/", contrast, "_", context, "_DMR_merged_summary.csv"))

####### 5. Retain context specific tile lists for comparison

CG_tiles_tested <- dmlTest_calls %>% mutate(tile = paste0(chr, "_", pos)) %>% select(tile)

CG_DMRs <- dmls_tbl %>% mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, 1, -1)) %>% select(tile, DMR)

CG_mu <- dmls_tbl %>%
  mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, "hyper", "hypo")) %>%
  mutate(DMR_type = paste0(context, ".", DMR)) %>%
  select(tile, DMR_type, mu1, mu2)

#######  #######
# CHG
#######  #######

context = "CHG"

#######  read in the files
sample1_data <- read_tsv(file.path(path_to_data_files, paste0(sample1, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample1_data <- sample1_data %>% select("chr", "pos", "N", "X", "sites")
# sample1_data <- sample1_data %>% filter(chr == "Chr01")
sample1_data

# # A tibble: 3,997,222 × 5
#    chr     pos     N     X sites
#    <chr> <dbl> <dbl> <dbl> <dbl>
#  1 Chr01   601    61    11     6
#  2 Chr01   701    64     0     5
#  3 Chr01   801   109     1     9
#  4 Chr01   901    87     0     7
#  5 Chr01  1301   100    49     5

sample2_data <- read_tsv(file.path(path_to_data_files, paste0(sample2, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample2_data <- sample2_data %>% select("chr", "pos", "N", "X", "sites")
# sample2_data <- sample2_data %>% filter(chr == "Chr01")
sample2_data

# # A tibble: 3,860,343 × 5
#    chr     pos     N     X sites
#    <chr> <dbl> <dbl> <dbl> <dbl>
#  1 Chr01   601    50     7     6
#  2 Chr01   701    53     0     5
#  3 Chr01   801   110     0     9
#  4 Chr01   901    66     0     7
#  5 Chr01  1301    63    27     5

#######  1. create an object of BSseq class (requires bsseq Bioconductor package)

BSobj <- makeBSseqData(list(sample1_data, sample2_data), c(sample1, sample2) )
BSobj

# An object of type 'BSseq' with
#   4044301 methylation loci
#   2 samples
# has not been smoothed
# All assays are in-memory

####### 2 Perform statistical test for DML by calling DMLtest function
# t1 <- proc.time()
dmlTest <- DMLtest(BSobj, group1=sample1, group2=sample2, smoothing = F, equal.disp = T, ncores = 4)
# proc.time() - t1

# user  system elapsed
# 508.827   0.139 509.181

####### 3. Call DMRs

dmls <- callDML(dmlTest, p.threshold=0.01, delta = 0.1)
dmls_tbl <- as.tibble(dmls)
dmls_tbl

# # A tibble: 4,480 × 12
#    chr        pos    mu1    mu2   diff diff.se  stat    phi1    phi2     pval
#    <chr>    <int>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>   <dbl>   <dbl>    <dbl>
#  1 Chr01  7169901 0.0685 0.586  -0.517  0.0671 -7.71 0.00674 0.00674 1.26e-14
#  2 Chr06 44281701 0.0723 0.551  -0.478  0.0626 -7.64 0.00674 0.00674 2.20e-14
#  3 Chr07  7827701 0.0434 0.562  -0.518  0.0698 -7.43 0.00674 0.00674 1.12e-13
#  4 Chr02 29499001 0.0948 0.570  -0.475  0.0647 -7.34 0.00674 0.00674 2.10e-13
#  5 Chr09 42493401 0.509  0.0469  0.462  0.0633  7.29 0.00674 0.00674 2.99e-13

# merge adjacent DMRs
# note: the bed must be ordered - tidy genomics does not check!
dmls_bed <- dmls_tbl %>% mutate(chr = as.character(chr), start = as.double(pos), end = pos+100, type = ifelse(diff > 0, "hyper", "hypo")) %>%
  select(chr, start, end, diff, type)  %>%
  arrange(chr, start)

dmls_bed
# # A tibble: 16,163 × 5
#    chr    start    end   diff type 
#    <chr>  <dbl>  <dbl>  <dbl> <chr>
#  1 Chr01   5701   5801  0.408 hyper
#  2 Chr01  35901  36001 -0.440 hypo 
#  3 Chr01  62101  62201  0.344 hyper
#  4 Chr01  78401  78501 -0.354 hypo 
#  5 Chr01  87301  87401  0.330 hyper

dmls_bed_summary <- dmls_bed %>%
  group_by(type) %>%
  mutate(length = end - start) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(diff))
dmls_bed_summary

# # A tibble: 2 × 5
#   type      n mean_length max_length mean_diff
#   <chr> <int>       <dbl>      <dbl>     <dbl>
# 1 hyper  2660         100        100     0.331
# 2 hypo   1820         100        100    -0.340

dmls_bed_hyper <- filter(dmls_bed, type == "hyper")
dmls_bed_hypo <- filter(dmls_bed, type == "hypo")

dmls_bed_hyper <- genome_cluster(dmls_bed_hyper, by=c("chr", "start", "end"), max_distance = 1)
dmls_bed_hypo <- genome_cluster(dmls_bed_hypo, by=c("chr", "start", "end"), max_distance = 1)

dmls_bed_pre_merge <- bind_rows(dmls_bed_hyper, dmls_bed_hypo)
# A tibble: 4,480 × 6

# dmls_bed %>% select(chr, start, end, diff) %>% arrange(chr, start)
# dmls_bed %>% select(chr, start, end, diff, cluster_id) %>% arrange(cluster_id)

dmls_bed_merged <- dmls_bed_pre_merge %>% group_by(cluster_id, type) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), mean_diff = mean(diff), max_diff = max(diff), min_diff = min(diff)) %>%
  mutate(length = end - start) %>%
  ungroup()
dmls_bed_merged

# # A tibble: 4,421 × 9
#    cluster_id type  chr      start      end mean_diff max_diff min_diff length
#         <dbl> <chr> <chr>    <dbl>    <dbl>     <dbl>    <dbl>    <dbl>  <dbl>
#  1          0 hyper Chr01   437901   438001     0.312    0.312    0.312    100
#  2          0 hypo  Chr01    20301    20401    -0.297   -0.297   -0.297    100
#  3          1 hyper Chr01   541301   541401     0.383    0.383    0.383    100
#  4          1 hypo  Chr01   267101   267201    -0.278   -0.278   -0.278    100
#  5          2 hyper Chr01  8179701  8179801     0.340    0.340    0.340    100

# dmls_bed_merged %>% filter(length > 100)

dmls_bed_merge_summarised <- dmls_bed_merged %>%
  group_by(type) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(mean_diff))
dmls_bed_merge_summarised

# # A tibble: 2 × 5
#   type      n mean_length max_length mean_diff
#   <chr> <int>       <dbl>      <dbl>     <dbl>
# 1 hyper  2626        101.        400     0.331
# 2 hypo   1795        101.        300    -0.341

####### 4. Make results files and write out tables

dmlTest_calls <- as.tibble(dmlTest) %>% left_join(dmls_tbl, by = colnames(dmlTest))
dmlTest_calls

### write DMR calls (large) file - only retain mu1 and mu2 cut off their columns to make file smaller
dmlTest_calls_out <- dmlTest_calls %>% select(chr, pos, mu1, mu2)
write_tsv(dmlTest_calls_out, paste0(outFolder, "/", contrast, "_", context, "_dmlTest_calls.tsv"))

### write only DMR tiles
write_csv(dmls_tbl, paste0(outFolder, "/", contrast, "_", context, "_DMRs.csv"))

### write summary of DMR tiles
write_csv(dmls_bed_summary, paste0(outFolder, "/", contrast, "_", context, "_DMR_tiles_summary.csv"))

### write merged DMRs
dmls_bed_merged_out <- dmls_bed_merged %>% select(chr, start, end, type, mean_diff, max_diff, min_diff, length)
write_tsv(dmls_bed_merged_out, paste0(outFolder, "/", contrast, "_", context, "_DMRs_merged.tsv"))

### write summary of merged DMR tiles
write_csv(dmls_bed_merge_summarised, paste0(outFolder, "/", contrast, "_", context, "_DMR_merged_summary.csv"))

####### 5. Retain contect sprcific tile lists for comparison

CHG_tiles_tested <- dmlTest_calls %>% mutate(tile = paste0(chr, "_", pos)) %>% select(tile)

CHG_DMRs <- dmls_tbl %>% mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, 1, -1)) %>% select(tile, DMR)

CHG_mu <- dmls_tbl %>%
  mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, "hyper", "hypo")) %>%
  mutate(DMR_type = paste0(context, ".", DMR)) %>%
  select(tile, DMR_type, mu1, mu2)

#######  #######
# CHH
#######  #######

context = "CHH"

#######  read in the files
sample1_data <- read_tsv(file.path(path_to_data_files, paste0(sample1, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample1_data <- sample1_data %>% select("chr", "pos", "N", "X", "sites")
# sample1_data <- sample1_data %>% filter(chr == "Chr01")
sample1_data

# # A tibble: 5,322,209 × 5
#    chr     pos     N     X sites
#    <chr> <dbl> <dbl> <dbl> <dbl>
#  1 Chr01   201   108     0    12
#  2 Chr01   301   393    25    28
#  3 Chr01   501   229     0    34
#  4 Chr01   601   195    19    19
#  5 Chr01   701   338     2    27

sample2_data <- read_tsv(file.path(path_to_data_files, paste0(sample2, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample2_data <- sample2_data %>% select("chr", "pos", "N", "X", "sites")
# sample2_data <- sample2_data %>% filter(chr == "Chr01")
sample2_data

# # A tibble: 5,098,155 × 5
#    chr     pos     N     X sites
#    <chr> <dbl> <dbl> <dbl> <dbl>
#  1 Chr01   201    65     0    12
#  2 Chr01   301   213     1    28
#  3 Chr01   601   186    24    19
#  4 Chr01   701   316     1    27
#  5 Chr01   801   368     1    32

#######  1. create an object of BSseq class (requires bsseq Bioconductor package)

BSobj <- makeBSseqData(list(sample1_data, sample2_data), c(sample1, sample2) )
BSobj

# An object of type 'BSseq' with
#   5381922 methylation loci
#   2 samples
# has not been smoothed
# All assays are in-memory

####### 2 Perform statistical test for DML by calling DMLtest function
# t1 <- proc.time()
dmlTest <- DMLtest(BSobj, group1=sample1, group2=sample2, smoothing = F, equal.disp = T, ncores = 4)
# proc.time() - t1

# user  system elapsed
# 508.827   0.139 509.181

####### 3. Call DMRs

dmls <- callDML(dmlTest, p.threshold=0.01, delta = 0.1)
dmls_tbl <- as.tibble(dmls)
dmls_tbl

# # A tibble: 1,171 × 12
#    chr        pos   mu1    mu2   diff diff.se  stat    phi1    phi2     pval
#    <chr>    <int> <dbl>  <dbl>  <dbl>   <dbl> <dbl>   <dbl>   <dbl>    <dbl>
#  1 Chr06 61115501 0.505 0.0307  0.474  0.0605  7.83 0.00674 0.00674 4.81e-15
#  2 Chr04 54573301 0.939 0.486   0.453  0.0583  7.78 0.00674 0.00674 7.41e-15
#  3 Chr04  7843001 0.557 0.0991  0.458  0.0598  7.66 0.00674 0.00674 1.81e-14
#  4 Chr04 57767601 0.690 0.217   0.472  0.0627  7.54 0.00674 0.00674 4.86e-14
#  5 Chr08 58053401 0.636 0.178   0.458  0.0609  7.52 0.00674 0.00674 5.60e-14

# merge adjacent DMRs
# note: the bed must be ordered - tidy genomics does not check!
dmls_bed <- dmls_tbl %>% mutate(chr = as.character(chr), start = as.double(pos), end = pos+100, type = ifelse(diff > 0, "hyper", "hypo")) %>%
  select(chr, start, end, diff, type)  %>%
  arrange(chr, start)

dmls_bed
# # A tibble: 1,171 × 5
#    chr     start     end  diff type 
#    <chr>   <dbl>   <dbl> <dbl> <chr>
#  1 Chr01    1101    1201 0.198 hyper
#  2 Chr01    1201    1301 0.231 hyper
#  3 Chr01  585901  586001 0.272 hyper
#  4 Chr01  850301  850401 0.229 hyper
#  5 Chr01 1156601 1156701 0.305 hyper

dmls_bed_summary <- dmls_bed %>%
  group_by(type) %>%
  mutate(length = end - start) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(diff))
dmls_bed_summary

# # A tibble: 2 × 5
#   type      n mean_length max_length mean_diff
#   <chr> <int>       <dbl>      <dbl>     <dbl>
# 1 hyper  1040         100        100     0.271
# 2 hypo    131         100        100    -0.271

dmls_bed_hyper <- filter(dmls_bed, type == "hyper")
dmls_bed_hypo <- filter(dmls_bed, type == "hypo")

dmls_bed_hyper <- genome_cluster(dmls_bed_hyper, by=c("chr", "start", "end"), max_distance = 1)
dmls_bed_hypo <- genome_cluster(dmls_bed_hypo, by=c("chr", "start", "end"), max_distance = 1)

dmls_bed_pre_merge <- bind_rows(dmls_bed_hyper, dmls_bed_hypo)
# A tibble: 1,171 × 6

# dmls_bed %>% select(chr, start, end, diff) %>% arrange(chr, start)
# dmls_bed %>% select(chr, start, end, diff, cluster_id) %>% arrange(cluster_id)

dmls_bed_merged <- dmls_bed_pre_merge %>% group_by(cluster_id, type) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), mean_diff = mean(diff), max_diff = max(diff), min_diff = min(diff)) %>%
  mutate(length = end - start) %>%
  ungroup()
dmls_bed_merged

# # A tibble: 1,140 × 9
#    cluster_id type  chr      start      end mean_diff max_diff min_diff length
#         <dbl> <chr> <chr>    <dbl>    <dbl>     <dbl>    <dbl>    <dbl>  <dbl>
#  1          0 hyper Chr01     1101     1301     0.215    0.231    0.198    200
#  2          0 hypo  Chr01  3850201  3850301    -0.278   -0.278   -0.278    100
#  3          1 hyper Chr01   585901   586001     0.272    0.272    0.272    100
#  4          1 hypo  Chr01  5063501  5063601    -0.247   -0.247   -0.247    100
#  5          2 hyper Chr01  2882501  2882601     0.326    0.326    0.326    10

# dmls_bed_merged %>% filter(length > 100)

dmls_bed_merge_summarised <- dmls_bed_merged %>%
  group_by(type) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(mean_diff))
dmls_bed_merge_summarised

# # A tibble: 2 × 5
#   type      n mean_length max_length mean_diff
#   <chr> <int>       <dbl>      <dbl>     <dbl>
# 1 hyper  1010        103.        200     0.271
# 2 hypo    130        101.        200    -0.271

####### 4. Make results files and write out tables

dmlTest_calls <- as.tibble(dmlTest) %>% left_join(dmls_tbl, by = colnames(dmlTest))
dmlTest_calls

### write DMR calls (large) file - only retain mu1 and mu2 cut off their columns to make file smaller
dmlTest_calls_out <- dmlTest_calls %>% select(chr, pos, mu1, mu2)
write_tsv(dmlTest_calls_out, paste0(outFolder, "/", contrast, "_", context, "_dmlTest_calls.tsv"))

### write only DMR tiles
write_csv(dmls_tbl, paste0(outFolder, "/", contrast, "_", context, "_DMRs.csv"))

### write summary of DMR tiles
write_csv(dmls_bed_summary, paste0(outFolder, "/", contrast, "_", context, "_DMR_tiles_summary.csv"))

### write merged DMRs
dmls_bed_merged_out <- dmls_bed_merged %>% select(chr, start, end, type, mean_diff, max_diff, min_diff, length)
write_tsv(dmls_bed_merged_out, paste0(outFolder, "/", contrast, "_", context, "_DMRs_merged.tsv"))

### write summary of merged DMR tiles
write_csv(dmls_bed_merge_summarised, paste0(outFolder, "/", contrast, "_", context, "_DMR_merged_summary.csv"))

####### 5. Retain contect sprcific tile lists for comparison

CHH_tiles_tested <- dmlTest_calls %>% mutate(tile = paste0(chr, "_", pos)) %>% select(tile)

CHH_DMRs <- dmls_tbl %>% mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, 1, -1)) %>% select(tile, DMR)

CHH_mu <- dmls_tbl %>%
  mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, "hyper", "hypo")) %>%
  mutate(DMR_type = paste0(context, ".", DMR)) %>%
  select(tile, DMR_type, mu1, mu2)

#######  ####### #######  ####### #######  ####### #######  ####### #######  ####### #######  #######
# Some analysis
#######  ####### #######  ####### #######  ####### #######  ####### #######  ####### #######  #######

outFolder_analysis <- paste0(outFolder, "/Analysis_", contrast)
dir.create(outFolder_analysis)

# Overlap between tiles that passed filtering

CG_tiles_tested <- CG_tiles_tested %>% mutate(context = "CG", pass = 1)
CHG_tiles_tested <- CHG_tiles_tested %>% mutate(context = "CHG", pass = 1)
CHH_tiles_tested <- CHH_tiles_tested %>% mutate(context = "CHH", pass = 1)

tile_lists <- bind_rows(CG_tiles_tested, CHG_tiles_tested, CHH_tiles_tested)

table_venn <- tile_lists %>%
  mutate(pass = as.integer(pass)) %>%
  spread(key = context, value = pass, fill = 0)

table_venn_summary <- table_venn %>%
  group_by(CG, CHG, CHH) %>%
  summarise(Counts = n()) %>%
  ungroup() %>%
  mutate(tiles_pct = round(Counts/sum(Counts)*100, digits = 1))

table_venn_summary

ref_venn <- as.tibble(expand.grid(0:1, 0:1, 0:1)) %>% mutate(Counts = 0, tiles_pct = 0)
colnames(ref_venn) <- colnames(table_venn_summary)

ref_venn

table_venn_summary_all <- ref_venn %>% anti_join(table_venn_summary, by = c(colnames(table_venn_summary)[1:3])) %>% bind_rows(table_venn_summary)

table_venn_summary_counts <- table_venn_summary_all %>%
  select(-tiles_pct)

venn_data <- as.matrix(table_venn_summary_counts)
class(venn_data) <- "VennCounts"
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_tiles_assessed_context_overlap.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "Tiles assessed per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

table_venn_summary_2 <- table_venn_summary_all %>%
  mutate(Counts = tiles_pct) %>%
  select(-tiles_pct)

venn_data <- as.matrix(table_venn_summary_2)
class(venn_data) <- "VennCounts"
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_tiles_assessed_context_overlap_pct.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "Tiles assessed per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

# Overlap between DMRs

CG_DMRs <- CG_DMRs %>% mutate(context = "CG")
CHG_DMRs <- CHG_DMRs %>% mutate(context = "CHG")
CHH_DMRs <- CHH_DMRs %>% mutate(context = "CHH")

tile_lists <- bind_rows(CG_DMRs, CHG_DMRs, CHH_DMRs)
tile_lists

table_venn <- tile_lists %>%
  mutate(DMR = as.integer(DMR)) %>%
  spread(key = context, value = DMR, fill = 0)

venn_data <- vennCounts(table_venn[,-1], include=c("both"))
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_both.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "all DMRs per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)

vennDiagram(venn_data, names = "",
            # main= "hyper DMRs per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()


venn_data <- vennCounts(table_venn[,-1], include=c("up"))
venn_data
pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_up.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "hyper DMRs per context",
            include=c("up"),
            counts.col=c('red'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)

vennDiagram(venn_data, names = "",
            # main= "hyper DMRs per context",
            include=c("up"),
            counts.col=c('red'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

venn_data <- vennCounts(table_venn[,-1], include=c("down"))
venn_data
pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_down.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "hypo DMRs per context",
            include=c("down"),
            counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)

vennDiagram(venn_data, names = "",
            # main= "hypo DMRs per context",
            include=c("down"),
            counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

venn_data <- vennCounts(table_venn[,-1], include=c("both"))
venn_data

venn_data_2 <- venn_data
class(venn_data_2) <- "matrix"

venn_data_pct <- as.tibble(venn_data_2) %>%
  mutate(Counts = round(Counts/sum(Counts)*100, digits = 1))

venn_data_pct

venn_data <- as.matrix(venn_data_pct)
class(venn_data) <- "VennCounts"
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_both_pct.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "DMRs per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

###### Density plot of mC levels for DMRs

DMR_lists <- bind_rows(CG_mu, CHG_mu, CHH_mu) %>% gather(key = sample, value = mC, -DMR_type, -tile)
DMR_lists

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_mC_values_distribution.pdf"), width = 6, height = 4)
ggplot(DMR_lists, aes(mC, colour = DMR_type)) +
  geom_density() +
  theme_minimal() +
  labs(title = "DMRs mC values distribution") +
  scale_color_ptol()
dev.off()

DMR_lists_2 <- DMR_lists %>% separate(DMR_type, into = c("context", "direction"))

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_mC_values_distribution_facet.pdf"), width = 10, height = 7)
ggplot(DMR_lists_2, aes(mC, colour = sample)) +
  geom_density() +
  theme_minimal() +
  labs(title = "DMRs mC values distribution") +
  scale_color_ptol() +
  facet_grid(direction ~ context)
dev.off()
