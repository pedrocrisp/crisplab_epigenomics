#!/usr/bin/Rscript
##########

#Peter Crisp
#2020-11-30
#R script to convert coverage matric to per bin average

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample_ID <- args[1]
sample_ID

out_dir <- args[1]
out_dir

# sample_ID = "UMR_McrBC_MF_merged"

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
library(ggplot2)
library(ggthemes)

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

###########################

#########  3. Merge in C sites data ##########

matrix_test <- read_tsv(paste0(out_dir, "/",sample_ID, ".mat_values.tab"), skip = 3, col_names = as.character(c(1:430)))

# replace NaN with 0
matrix_test_zero <- matrix_test %>%
  # slice(1:10000) %>%
  # select(-UMR_McrBC_MF_merged_429) %>%
  replace(., is.na(.), 0)

# gather and summarise
matrix_test_zero_gather <- matrix_test_zero %>%
  gather(key = bin, value = coverage) %>%
  mutate(bin = as.double(bin)) %>%
  group_by(bin) %>%
  summarise(average_coverage = mean(coverage))

matrix_test_zero_gather

write.table(reference_tiles_sites, paste0(out_dir, "/", sample_ID, "mat_values_summarised.tab"), sep = "\t", quote = F, row.names = F)


g <- ggplot(matrix_test_zero_gather, aes(x = bin, y = average_coverage)) +
  geom_line() +
  theme_classic() +
  text_size_theme_8
g

ggsave(filename =paste0(out_dir, "/", sample_ID, "mat_metaplot_2.pdf"), plot =  g, width = 3, height = 3)
