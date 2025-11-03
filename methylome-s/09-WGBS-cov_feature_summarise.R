#!/usr/bin/Rscript
##########

#Peter Crisp
#2021-06-21
#R script to convert coverage matrix to per bin average for each of CG, CHG and CHH

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample_ID <- args[1]
sample_ID

out_dir <- args[2]
out_dir

# sample_ID = "maize_0.5ug_120m_S29"
# out_dir = "trimmed_align_bowtie2_bigWigs_deeptools_WGBS-levels-100-bp"

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

######### CG ########
context = "CG"
matrix_test <- read_tsv(paste0(out_dir, "/",sample_ID, "_", context,".mat_values.tab"), skip = 3, col_names = as.character(c(1:43)))

# DONT replace NaN with 0 for mC data!
matrix_test_zero <- matrix_test
  # %>%
  # slice(1:10000) %>%
  # select(-UMR_McrBC_MF_merged_429) %>%
  # replace(., is.na(.), 0)

# gather and summarise
matrix_test_zero_gather <- matrix_test_zero %>%
  gather(key = bin, value = coverage) %>%
  mutate(bin = as.double(bin)) %>%
  group_by(bin) %>%
  summarise(average_coverage = mean(coverage, na.rm = T))

matrix_test_zero_gather

write.table(matrix_test_zero_gather, paste0(out_dir, "/", sample_ID, "_", context,".mat_values_summarised.tab"), sep = "\t", quote = F, row.names = F)


g <- ggplot(matrix_test_zero_gather, aes(x = bin, y = average_coverage)) +
  geom_line() +
  theme_classic() +
  text_size_theme_8
g

ggsave(filename =paste0(out_dir, "/", sample_ID, "_", context, ".mat_metaplot_2.pdf"), plot =  g, width = 3, height = 3)

######### CHG ########
context = "CHG"
matrix_test <- read_tsv(paste0(out_dir, "/",sample_ID, "_", context,".mat_values.tab"), skip = 3, col_names = as.character(c(1:43)))

# DONT replace NaN with 0 for mC data!
matrix_test_zero <- matrix_test
# %>%
# slice(1:10000) %>%
# select(-UMR_McrBC_MF_merged_429) %>%
# replace(., is.na(.), 0)

# gather and summarise
matrix_test_zero_gather <- matrix_test_zero %>%
  gather(key = bin, value = coverage) %>%
  mutate(bin = as.double(bin)) %>%
  group_by(bin) %>%
  summarise(average_coverage = mean(coverage, na.rm = T))

matrix_test_zero_gather

write.table(matrix_test_zero_gather, paste0(out_dir, "/", sample_ID, "_", context,".mat_values_summarised.tab"), sep = "\t", quote = F, row.names = F)


g <- ggplot(matrix_test_zero_gather, aes(x = bin, y = average_coverage)) +
  geom_line() +
  theme_classic() +
  text_size_theme_8
g

ggsave(filename =paste0(out_dir, "/", sample_ID, "_", context, ".mat_metaplot_2.pdf"), plot =  g, width = 3, height = 3)

######### CHH ########
context = "CHH"
matrix_test <- read_tsv(paste0(out_dir, "/",sample_ID, "_", context,".mat_values.tab"), skip = 3, col_names = as.character(c(1:43)))

# DONT replace NaN with 0 for mC data!
matrix_test_zero <- matrix_test
# %>%
# slice(1:10000) %>%
# select(-UMR_McrBC_MF_merged_429) %>%
# replace(., is.na(.), 0)

# gather and summarise
matrix_test_zero_gather <- matrix_test_zero %>%
  gather(key = bin, value = coverage) %>%
  mutate(bin = as.double(bin)) %>%
  group_by(bin) %>%
  summarise(average_coverage = mean(coverage, na.rm = T))

matrix_test_zero_gather

write.table(matrix_test_zero_gather, paste0(out_dir, "/", sample_ID, "_", context,".mat_values_summarised.tab"), sep = "\t", quote = F, row.names = F)


g <- ggplot(matrix_test_zero_gather, aes(x = bin, y = average_coverage)) +
  geom_line() +
  theme_classic() +
  text_size_theme_8
g

ggsave(filename =paste0(out_dir, "/", sample_ID, "_", context, ".mat_metaplot_2.pdf"), plot =  g, width = 3, height = 3)

