#!/usr/bin/Rscript
##########

#Peter Crisp
#2023-10-10
#R script to determine combine mC context to single vlaue per tiles, 
# using the max mC of the contexts per tile
# suggest using a filter of 2 sites and coverage 3

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample_name <- args[1]
outFolder <- args[2]
coverage_filter_min <- as.double(args[3])
site_filter_min <- as.double(args[4])

#args
# sample_name = "maize_B73_leaf_EM-2-3"
# outFolder = "tiles_C"

########
###########################
library(tidyverse)
old.scipen <- getOption("scipen")
options(scipen=999)

dir.create(outFolder, showWarnings = F)

cg <- read_tsv(paste0("tiles/", "sample_name", "_merged_BSMAP_out.txt.100.CG.fixed.sorted.txt"))
chg <- read_tsv(paste0("tiles/", "sample_name", "_merged_BSMAP_out.txt.100.CHG.fixed.sorted.txt"))
chh <- read_tsv(paste0("tiles/", "sample_name", "_merged_BSMAP_out.txt.100.CHH.fixed.sorted.txt"))

cg_filtered <- cg %>% 
  mutate(cov = CT/cg_sites) %>% 
  filter(cg_sites >= site_filter_min & cov >= coverage_filter_min) %>%
  select(chr, start, end, ratio) %>%
  mutate(context = "CG")
# A tibble: 15,230,273 × 9
chg_filtered <- chg %>% 
  mutate(cov = CT/chg_sites) %>% 
  filter(chg_sites >= site_filter_min & cov >= coverage_filter_min) %>%
  select(chr, start, end, ratio) %>%
  mutate(context = "CHG")
# A tibble: 15,569,899 × 9
chh_filtered <- chh %>% 
  mutate(cov = CT/chh_sites) %>% 
  filter(chh_sites >= site_filter_min & cov >= coverage_filter_min) %>%
  select(chr, start, end, ratio) %>%
  mutate(context = "CHH")
# A tibble: 17,529,863 × 9

c <- bind_rows(cg_filtered, chg_filtered, chh_filtered)

# # A tibble: 48,330,035 × 5
#    chr   start   end ratio context
#    <chr> <dbl> <dbl> <dbl> <chr>  
#  1 1       101   200 0.825 CG     
#  2 1       201   300 1     CG     
#  3 1       301   400 0.857 CG    

c_max <- c %>% group_by(chr, start, end) %>%
  summarise(ratio = max(ratio)) %>% 
  ungroup() %>% mutate(start = start-1) 

# # A tibble: 17,980,874 × 4
#    chr   start   end  ratio
#    <chr> <dbl> <dbl>  <dbl>
#  1 1         0   100 0.0178
#  2 1       100   200 0.825 
#  3 1       200   300 1     
#  4 1       300   400 0.857 

write.table(c_max, paste0(outFolder,"/", "sample_name","_merged_BSMAP_out.txt.100.C.max.bed.txt"), sep = "\t", quote = F, row.names = F, col.names = T)