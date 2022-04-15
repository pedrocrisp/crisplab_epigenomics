#!/usr/bin/Rscript
##########

# Peter Crisp
# 2022-04-03
# R script to summarise MethylDackel output into 100bp tiles
# then
# R script to amend 100bp tile data with incorrect chr end tiles
# because 100bp runs past the end of the chr in most cases

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample <- args[1]
sample
data_folder <- args[2]
data_folder
reference_tile_file <- args[3]
reference_tile_file

# ##### dev
# cd /scratch/project/crisp006/pete/arabidopsis-bwa-WGBS/analysis
# module load R/3.5.0-gnu
# mkdir -p tiles
# R
# sample <- "SRR5724529_Chr1"
# # SRR5724529_Chr1_methratio_CG.bedGraph
# data_folder <- "MethylDackel"
# reference_tile_file <- "/home/uqpcrisp/refseqs/arabidopsis/chromosomes/sites/Chr1_100bp_tiles_zBased_sites_counts.txt"

###########################
#setup
library(tidyverse)
library(RcppRoll)
library(purrr)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)

outdir = "tiles"
###########################

#reference used to make amendments
reference_tiles <- read_tsv(reference_tile_file, col_names = TRUE,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer(),
                              start_zBased = col_integer(),
                              cg_sites = col_integer(),
                              chg_sites = col_integer(),
                              chh_sites = col_integer()
                            ))


### ### ### ### ### ###
# Make tiles module CG
### ### ### ### ### ###

mc_tiles <- read_tsv(paste0(data_folder, "/", sample, "_methratio_CG.bedGraph"),
                     col_names = c("chr", "start", "end", "ratio", "C", "CT"))
mc_tiles

# # A tibble: 1,415,283 × 6
# chr   start   end ratio     C    CT
# <chr> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Chr1     31    32   100     1     0
# 2 Chr1     33    34    50     5     5
# 3 Chr1    107   108    50     5     5
# 4 Chr1    113   114    22     2     7
# 5 Chr1    116   117    20     3    12
# 6 Chr1    122   123    66     6     3
# 7 Chr1    124   125    86    13     2
# 8 Chr1    388   389     0     0     7
# 9 Chr1    390   391     7     1    13
# 10 Chr1    561   562     0     0     9

#"sites_with_data", "C", "CT", "ratio"

broken_bedGraph <- mc_tiles %>%
  mutate(start_zBased = start - start %% 100,
         broken_end = start_zBased + 100) %>%
  group_by(chr, start_zBased, broken_end) %>%
  summarise(sites_with_data = n(),
            C = sum(C),
            CT = sum(CT),
            ratio = mean(ratio)) %>% # summarise to the bins
  arrange(chr, start_zBased)

broken_bedGraph

# # A tibble: 256,629 × 7
# # Groups:   chr, start_zBased [256,629]
# chr   start_zBased broken_end sites_with_data     C    CT ratio
# <chr>        <dbl>      <dbl>           <int> <dbl> <dbl> <dbl>
# 1 Chr1             0        100               2     6     5 75
# 2 Chr1           100        200               5    29    29 48.8
# 3 Chr1           300        400               2     1    20  3.5
# 4 Chr1           500        600               4     3    37  7.5
# 5 Chr1           700        800               3    36    11 77.7
# 6 Chr1           800        900               5    55    26 66
# 7 Chr1           900       1000               3    16    23 43
# 8 Chr1          1000       1100               2     0     6  0
# 9 Chr1          1100       1200               2     0    27  0
# 10 Chr1          1200       1300               4     1    58  1.25

# write.table(mc_tiles_mean, paste0("tiles_tmp", "/", sample, "_methratio_CG.100.txt"), sep = "\t", quote = F, row.names = F, col.names = T)


### ### ### ### ### ### ###
# ADMEND TILE END MODULE CG
### ### ### ### ### ### ###

#########
# fix bed CG
# note: using write.table, write_delim converts to scientific notation

# # and cannot seem to disable that
# broken_bedGraph <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.bed"), col_names = F,
#                             cols(
#                               X1 = col_character(),
#                               X2 = col_integer(),
#                               X3 = col_integer(),
#                               X4 = col_number()
#                             ))
#
# colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "sites_with_data", "C", "CT", "ratio")
# #remove NAs - this is necessary because the output from the perl script is any tile with data in any context per sample
# broken_bedGraph <- na.omit(broken_bedGraph)

#fix
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))

# # uncomment if you want to make a bigwig file later
# #keep as BED format (zero-based coordinates)
# # only retain chr, start, end, ratio; then sort to make it a bedGraph file
# fixed_bedGraph2 <-
# fixed_bedGraph %>% select(chr, start_zBased, end, ratio) %>% arrange(chr, start_zBased)
# write.table(fixed_bedGraph2, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)
# rm(fixed_bedGraph2)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph3 <-
fixed_bedGraph %>% select(chr, start, end, C, CT, ratio, sites_with_data, cg_sites) %>% arrange(chr, start)
head(fixed_bedGraph3)
write.table(fixed_bedGraph3, paste0(outdir, "/", sample, "_BSMAP_out.txt.100.CG.fixed.sorted.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

# clear memory
rm(fixed_bedGraph3)
rm(fixed_bedGraph)
rm(broken_bedGraph)

### ### ### ### ### ###
# Make tiles module CHG
### ### ### ### ### ###

mc_tiles <- read_tsv(paste0(data_folder, "/", sample, "_methratio_CHG.bedGraph"),
                     col_names = c("chr", "start", "end", "ratio", "C", "CT"))
mc_tiles

# # A tibble: 1,415,283 × 6
# chr   start   end ratio     C    CT
# <chr> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Chr1     31    32   100     1     0
# 2 Chr1     33    34    50     5     5
# 3 Chr1    107   108    50     5     5
# 4 Chr1    113   114    22     2     7
# 5 Chr1    116   117    20     3    12
# 6 Chr1    122   123    66     6     3
# 7 Chr1    124   125    86    13     2
# 8 Chr1    388   389     0     0     7
# 9 Chr1    390   391     7     1    13
# 10 Chr1    561   562     0     0     9

#"sites_with_data", "C", "CT", "ratio"

broken_bedGraph <- mc_tiles %>%
  mutate(start_zBased = start - start %% 100,
         broken_end = start_zBased + 100) %>%
  group_by(chr, start_zBased, broken_end) %>%
  summarise(sites_with_data = n(),
            C = sum(C),
            CT = sum(CT),
            ratio = mean(ratio)) %>% # summarise to the bins
  arrange(chr, start_zBased)

broken_bedGraph

# # A tibble: 256,629 × 7
# # Groups:   chr, start_zBased [256,629]
# chr   start_zBased broken_end sites_with_data     C    CT ratio
# <chr>        <dbl>      <dbl>           <int> <dbl> <dbl> <dbl>
# 1 Chr1             0        100               2     6     5 75
# 2 Chr1           100        200               5    29    29 48.8
# 3 Chr1           300        400               2     1    20  3.5
# 4 Chr1           500        600               4     3    37  7.5
# 5 Chr1           700        800               3    36    11 77.7
# 6 Chr1           800        900               5    55    26 66
# 7 Chr1           900       1000               3    16    23 43
# 8 Chr1          1000       1100               2     0     6  0
# 9 Chr1          1100       1200               2     0    27  0
# 10 Chr1          1200       1300               4     1    58  1.25

# write.table(mc_tiles_mean, paste0("tiles_tmp", "/", sample, "_methratio_CHG.100.txt"), sep = "\t", quote = F, row.names = F, col.names = T)


### ### ### ### ### ### ###
# ADMEND TILE END MODULE CHG
### ### ### ### ### ### ###

#########
# fix bed CHG
# note: using write.table, write_delim converts to scientific notation

# # and cannot seem to disable that
# broken_bedGraph <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.bed"), col_names = F,
#                             cols(
#                               X1 = col_character(),
#                               X2 = col_integer(),
#                               X3 = col_integer(),
#                               X4 = col_number()
#                             ))
#
# colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "sites_with_data", "C", "CT", "ratio")
# #remove NAs - this is necessary because the output from the perl script is any tile with data in any context per sample
# broken_bedGraph <- na.omit(broken_bedGraph)

#fix
fixed_bedGraph <-
  merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))

# # uncomment if you want to make a bigwig file later
# #keep as BED format (zero-based coordinates)
# # only retain chr, start, end, ratio; then sort to make it a bedGraph file
# fixed_bedGraph2 <-
# fixed_bedGraph %>% select(chr, start_zBased, end, ratio) %>% arrange(chr, start_zBased)
# write.table(fixed_bedGraph2, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)
# rm(fixed_bedGraph2)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph3 <-
  fixed_bedGraph %>% select(chr, start, end, C, CT, ratio, sites_with_data, chg_sites) %>% arrange(chr, start)
write.table(fixed_bedGraph3, paste0(outdir, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.sorted.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

# clear memory
rm(fixed_bedGraph3)
rm(fixed_bedGraph)
rm(broken_bedGraph)


### ### ### ### ### ###
# Make tiles module CHH
### ### ### ### ### ###

mc_tiles <- read_tsv(paste0(data_folder, "/", sample, "_methratio_CHH.bedGraph"),
                     col_names = c("chr", "start", "end", "ratio", "C", "CT"))
mc_tiles

# # A tibble: 1,415,283 × 6
# chr   start   end ratio     C    CT
# <chr> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Chr1     31    32   100     1     0
# 2 Chr1     33    34    50     5     5
# 3 Chr1    107   108    50     5     5
# 4 Chr1    113   114    22     2     7
# 5 Chr1    116   117    20     3    12
# 6 Chr1    122   123    66     6     3
# 7 Chr1    124   125    86    13     2
# 8 Chr1    388   389     0     0     7
# 9 Chr1    390   391     7     1    13
# 10 Chr1    561   562     0     0     9

#"sites_with_data", "C", "CT", "ratio"

broken_bedGraph <- mc_tiles %>%
  mutate(start_zBased = start - start %% 100,
         broken_end = start_zBased + 100) %>%
  group_by(chr, start_zBased, broken_end) %>%
  summarise(sites_with_data = n(),
            C = sum(C),
            CT = sum(CT),
            ratio = mean(ratio)) %>% # summarise to the bins
  arrange(chr, start_zBased)

broken_bedGraph

# # A tibble: 256,629 × 7
# # Groups:   chr, start_zBased [256,629]
# chr   start_zBased broken_end sites_with_data     C    CT ratio
# <chr>        <dbl>      <dbl>           <int> <dbl> <dbl> <dbl>
# 1 Chr1             0        100               2     6     5 75
# 2 Chr1           100        200               5    29    29 48.8
# 3 Chr1           300        400               2     1    20  3.5
# 4 Chr1           500        600               4     3    37  7.5
# 5 Chr1           700        800               3    36    11 77.7
# 6 Chr1           800        900               5    55    26 66
# 7 Chr1           900       1000               3    16    23 43
# 8 Chr1          1000       1100               2     0     6  0
# 9 Chr1          1100       1200               2     0    27  0
# 10 Chr1          1200       1300               4     1    58  1.25

# write.table(mc_tiles_mean, paste0("tiles_tmp", "/", sample, "_methratio_CHH.100.txt"), sep = "\t", quote = F, row.names = F, col.names = T)


### ### ### ### ### ### ###
# ADMEND TILE END MODULE CHH
### ### ### ### ### ### ###

#########
# fix bed CHH
# note: using write.table, write_delim converts to scientific notation

# # and cannot seem to disable that
# broken_bedGraph <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.bed"), col_names = F,
#                             cols(
#                               X1 = col_character(),
#                               X2 = col_integer(),
#                               X3 = col_integer(),
#                               X4 = col_number()
#                             ))
#
# colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "sites_with_data", "C", "CT", "ratio")
# #remove NAs - this is necessary because the output from the perl script is any tile with data in any context per sample
# broken_bedGraph <- na.omit(broken_bedGraph)

#fix
fixed_bedGraph <-
  merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))

# # uncomment if you want to make a bigwig file later
# #keep as BED format (zero-based coordinates)
# # only retain chr, start, end, ratio; then sort to make it a bedGraph file
# fixed_bedGraph2 <-
# fixed_bedGraph %>% select(chr, start_zBased, end, ratio) %>% arrange(chr, start_zBased)
# write.table(fixed_bedGraph2, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)
# rm(fixed_bedGraph2)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph3 <-
  fixed_bedGraph %>% select(chr, start, end, C, CT, ratio, sites_with_data, chh_sites) %>% arrange(chr, start)
write.table(fixed_bedGraph3, paste0(outdir, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.sorted.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

# clear memory
rm(fixed_bedGraph3)
rm(fixed_bedGraph)
rm(broken_bedGraph)
