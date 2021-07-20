#!/usr/bin/Rscript
##########

# Peter Crisp
# 2021-07-20
# R script to summary stats for 100bp mC tile data

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample <- args[1]
sample
data_folder <- args[2]
data_folder
out_folder <- args[3]
out_folder
out_folder <- args[4]
out_folder

###########################
#setup
library(tidyverse)
library(purrr)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
###########################

# args debugging
# sample = "Sb_mc_01"
# data_folder = "tiles_filtered_4C_5x"
# out_folder = "tiles_filtered_4C_5x_tileMeans"
# dir.create(showWarnings = F, out_folder)

###### CG ######
context = "CG"
tile_file = read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.", context, "_filtered.txt"))

tile_file_summary <- tile_file %>% mutate(ratio = C/CT*100) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

tile_file_summary

write.table(tile_file_summary, paste0(out_folder, "/", sample, "_", context, "_tile_means.tsv"), sep = "\t", quote = F, row.names = F)

###### CHG ######
context = "CHG"
tile_file = read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.", context, "_filtered.txt"))

tile_file_summary <- tile_file %>% mutate(ratio = C/CT*100) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

tile_file_summary

write.table(tile_file_summary, paste0(out_folder, "/", sample, "_", context, "_tile_means.tsv"), sep = "\t", quote = F, row.names = F)

###### CHH ######
context = "CHH"
tile_file = read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.", context, "_filtered.txt"))

tile_file_summary <- tile_file %>% mutate(ratio = C/CT*100) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

tile_file_summary

write.table(tile_file_summary, paste0(out_folder, "/", sample, "_", context, "_tile_means.tsv"), sep = "\t", quote = F, row.names = F)
