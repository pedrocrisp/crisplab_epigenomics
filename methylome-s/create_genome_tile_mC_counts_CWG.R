#!/usr/bin/Rscript
##########

#Peter Crisp
#2019-12-12
#R script to count number of CHG subcontext (CAG, CTG, CCG) sites per 100bp tile

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
genome <- args[1]
genome

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

# genome = "Sbicolor_454_v3.0.1"
reference_tiles <- read_tsv(paste0(genome, "_100bp_tiles_zBased.txt"), col_names = TRUE,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer(),
                              start_zBased = col_integer()
                            ))

unique(reference_tiles$chr)

reference_tiles_2 <- reference_tiles
# %>% filter(chr %in% c(1:10, "Pt", "Mt"))

reference_tiles_2
# A tibble: 7,097,206 x 4
#    chr   start   end start_zBased
#    <chr> <int> <int>        <int>
#  1 Chr11     1   100            0
#  2 Chr11   101   200          100
#  3 Chr11   201   300          200
#  4 Chr11   301   400          300
#  5 Chr11   401   500          400
#  6 Chr11   501   600          500
#  7 Chr11   601   700          600
#  8 Chr11   701   800          700
#  9 Chr11   801   900          800
# 10 Chr11   901  1000          900

# unique(reference_tiles_2$chr)

##### read in site files and parse

# read in CG/CHG/CHH sites reference file created above
reference_tiles_C_sites <- read_tsv(paste0(genome, "_100bp_tiles_sites.tsv"), col_names = TRUE,
                                    cols(
                                      seqID = col_character(),
                                      patternName = col_character(),
                                      frequency = col_double()
                                    ))

# parse to make a bed file
reference_tiles_C_sites_bed <- reference_tiles_C_sites %>%
  #head(n=1000) %>%
  separate(seqID, into = c("chr", "start", "end"), sep = c(":|-")) %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  mutate(start = ifelse(start >= 1, start+2, 0), end = end-2) %>%
  spread(key = patternName, value = frequency) %>%
  arrange(chr, start) %>%
  rename(cag_sites = CAG,
         ctg_sites = CTG,
         ccg_sites = CCG
  )

reference_tiles_C_sites_bed

# # A tibble: 6,253,575 x 6
#    chr   start   end cg_sites chg_sites chh_sites
#    <chr> <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#  1 Chr00     0   100       12         8        34
#  2 Chr00   100   200       NA        NA        18
#  3 Chr00   200   300        4         3         8
#  4 Chr00   300   400       10         2        26
#  5 Chr00   400   500       10         4        26
#  6 Chr00   500   600       12         2        27
#  7 Chr00   600   700        8         2        25
#  8 Chr00   700   800        4         4        23
#  9 Chr00   800   900       NA         4        26
# 10 Chr00   900  1000       NA         4        31

# unique(reference_tiles_C_sites_bed$chr)

#merge
reference_tiles_sites <- left_join(reference_tiles_2, reference_tiles_C_sites_bed, by = c("chr", "start_zBased" ="start", "end"))
reference_tiles_sites

# # A tibble: 7,097,206 x 7
#    chr   start   end start_zBased cg_sites chg_sites chh_sites
#    <chr> <int> <dbl>        <dbl>    <dbl>     <dbl>     <dbl>
#  1 Chr11     1   100            0        2         6        30
#  2 Chr11   101   200          100        6         3        29
#  3 Chr11   201   300          200        2         3        32
#  4 Chr11   301   400          300        4         3        33
#  5 Chr11   401   500          400        2         6        27
#  6 Chr11   501   600          500       NA         8        31
#  7 Chr11   601   700          600        2         5        21
#  8 Chr11   701   800          700        2         8        26
#  9 Chr11   801   900          800        2         6        25
# 10 Chr11   901  1000          900        2        11        37

# replace NA with 0 and re order
# which(is.na(reference_tiles_sites))

reference_tiles_sites <- reference_tiles_sites %>%
  replace(., is.na(.), 0) %>%
  mutate(chr = factor(chr, levels = unique(chr))) %>%
  arrange(chr, start)

unique(reference_tiles_sites$chr)

#write the one-based txt file
write.table(reference_tiles_sites, paste0(genome, "_100bp_tiles_zBased_sites_counts.txt"), sep = "\t", quote = F, row.names = F)

##############

#write a bed file (zero based and no head)
reference_tiles_2_bed <- select(reference_tiles_sites, chr, start_zBased, end, cag_sites, ctg_sites, cgg_sites)

write.table(reference_tiles_2_bed, paste0(genome, "_100bp_tiles_zBased_sites_counts.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
