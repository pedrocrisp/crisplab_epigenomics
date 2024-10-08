#!/usr/bin/Rscript
##########

#Peter Crisp
#2020-20-7
#R script to call UMRs using a 100bp tile

# Notes
# Currently this script removes the organelles by using filter(!chr %in% c("Mt", "Pt")) - this may not catch all organelles
# post filtering may be required to remove organelles or other undesired contigs

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
reference_100bp_tiles <- args[1]
sample_to_crunch <- args[2]
annotation_suffix <- args[3]
chrom_sizes_path <- args[4]
coverage_filter_min <- as.double(args[5])
site_filter_min <- as.double(args[6])
MR_percent <- as.double(args[7])
UMR_percent <- as.double(args[8])

######## de bug
# args
# reference_100bp_tiles = "/home/springer/pcrisp/ws/refseqs/maize/sites/Zea_mays.AGPv4_100bp_tiles_zBased_sites_counts.txt"
# sample_to_crunch = "generations_merged"
#
# # filter args
# coverage_filter_min = 3
# site_filter_min = 2
# MR_percent = 0.4
# UMR_percent = 0.1
# annotation_suffix = paste0("_mC_domains_II",
#                           "_cov_",coverage_filter_min,
#                           "_sites_",site_filter_min,
#                           "_MR_",MR_percent,
#                           "_UMR_",UMR_percent)
########
###########################
library(tidyverse)
library(ggthemes)
# library("seqinr")
old.scipen <- getOption("scipen")
options(scipen=999)
# library(wesanderson)
library(RColorBrewer)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
                           axis.title=element_text(size=8),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.title=element_text(size=8),
                           legend.text=element_text(size=8))

sample_to_crunch

###########################
# Module #1
###########################

folder_prefix = paste0("UMR_tiles_per_chr/", sample_to_crunch, annotation_suffix)

dir.create(folder_prefix, showWarnings = F)

out_dir = paste0(folder_prefix, "/mC_UMT_annotation")
dir.create(out_dir, showWarnings = F)
out_dir

out_dir_beds = paste0(folder_prefix, "/mC_UMT_annotation_beds")
dir.create(out_dir_beds, showWarnings = F)

######

#############
# folders
umr_out_dir = paste0(folder_prefix, "/mC_UMR_annotation")
dir.create(out_dir, showWarnings = F)

umr_out_dir_beds = paste0(folder_prefix, "/mC_UMR_annotation_beds")
dir.create(out_dir_beds, showWarnings = F)

#############
## read in mC tile data
#############

reference_tiles <- read_tsv(reference_100bp_tiles, col_names = T,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer()))

reference_tiles <- reference_tiles %>% select(-start_zBased)

reference_tiles %>% distinct(chr)

# Make decision about whether using contigs
reference_tiles %>% mutate(size = end - start) %>% summarise(MB = sum(size)/1000000)
# 702

###########

# read in CG
CG <- read_tsv(paste0("tiles/", sample_to_crunch, "_BSMAP_out.txt.100.CG.fixed.sorted.txt"), col_names = T, cols_only(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  CT = col_integer(),
#  C = col_integer(),
  ratio = col_double(),
#  cg_sites = col_integer())) %>% mutate(cov = (CT+C)/cg_sites)
  cg_sites = col_integer())) %>% mutate(cov = CT/cg_sites)

CG

g <- ggplot(CG, aes(x = cov)) +
  geom_density() +
  theme_minimal() +
  xlim(0,25) +
  geom_vline(xintercept = c(coverage_filter_min), colour = 'blue', linetype = 'dashed') +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/mC_CG_cov_density.pdf"), h = 1.5, w = 2)

# read in CHG
CHG <- read_tsv(paste0("tiles/", sample_to_crunch, "_BSMAP_out.txt.100.CHG.fixed.sorted.txt"), col_names = T, cols_only(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  CT = col_integer(),
#  C = col_integer(),
  ratio = col_double(),
#  chg_sites = col_integer())) %>% mutate(cov = (CT+C)/chg_sites)
  chg_sites = col_integer())) %>% mutate(cov = CT/chg_sites)

g <- ggplot(CHG, aes(x = cov)) +
  geom_density() +
  theme_minimal() +
  xlim(0,25) +
  geom_vline(xintercept = c(coverage_filter_min), colour = 'blue', linetype = 'dashed') +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/mC_CHG_cov_density.pdf"), h = 1.5, w = 2)

# read in CHH
CHH <- read_tsv(paste0("tiles/", sample_to_crunch, "_BSMAP_out.txt.100.CHH.fixed.sorted.txt"), col_names = T, cols_only(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  CT = col_integer(),
#  C = col_integer(),
  ratio = col_double(),
#  chh_sites = col_integer())) %>% mutate(cov = (CT+C)/chh_sites)
  chh_sites = col_integer())) %>% mutate(cov = CT/chh_sites)

g <- ggplot(CHH, aes(x = cov)) +
  geom_density() +
  theme_minimal() +
  xlim(0,25) +
  geom_vline(xintercept = c(coverage_filter_min), colour = 'blue', linetype = 'dashed') +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/mC_CHH_cov_density.pdf"), h = 1.5, w = 2)

########### ########### ###########
## Filter
########### ########### ###########

# filter
CG_ratio <- CG %>%
  mutate(CG = ifelse(cg_sites < site_filter_min | cov < coverage_filter_min, NA, ratio)) %>%
  select(chr, start, end, CG)

# filter
CHG_ratio <- CHG %>%
  mutate(CHG = ifelse(chg_sites < site_filter_min | cov < coverage_filter_min, NA, ratio)) %>%
  select(chr, start, end, CHG)

# filter
CHH_ratio <- CHH %>%
  mutate(CHH = ifelse(chh_sites < site_filter_min | cov < coverage_filter_min, NA, ratio)) %>%
  select(chr, start, end, CHH)

############
# merge with reference
# call tiles with no sites
merged_mC <- reference_tiles %>%
  # select("chr", "start", "end") %>%
  left_join(CG_ratio, by = c("chr", "start", "end")) %>%
  left_join(CHG_ratio, by = c("chr", "start", "end")) %>%
  left_join(CHH_ratio, by = c("chr", "start", "end")) %>%
  mutate(cg_sites = ifelse(cg_sites < site_filter_min, "n", "y"),
         chg_sites = ifelse(chg_sites < site_filter_min, "n", "y"),
         chh_sites = ifelse(chh_sites < site_filter_min, "n", "y"))

print(reference_tiles, n = 20)
print(CHG_ratio, n = 20)
print(merged_mC, n = 20)


# remove organelles
merged_mC %>% distinct(chr)
merged_mC %>% distinct(chr) %>% filter(grepl(paste("M", "P", sep = "|"), chr))
merged_mC_sans_orgs <- merged_mC %>% filter(!chr %in% c("Mt", "Pt"))

merged_mC_sans_orgs

###################
# distribution and averages of mC levels

n_for_sample_n = ifelse(nrow(merged_mC_sans_orgs) > 1000000, 1000000, nrow(merged_mC_sans_orgs))

# size of dataframe
n_for_sample_n

# distro
plot_data <- merged_mC_sans_orgs %>% 
sample_n(n_for_sample_n) %>%
  select(CG:CHH) %>%
  gather(key = context, value = percent) %>%
  mutate(percent = percent *100)

g <- ggplot(plot_data, aes(x = percent)) +
  geom_density() +
  geom_vline(xintercept = c(UMR_percent*100, MR_percent*100), colour = 'blue', linetype = 'dashed') +
  facet_grid(context ~., scales = 'free') +
  theme_minimal() +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/mC_tile_density.pdf"), h = 4, w = 2)

# mean
plot_data_average <- plot_data %>%
  group_by(context) %>%
  summarise(mean = mean(percent, na.rm = T))

g <- ggplot(plot_data_average, aes(x = context, y = mean)) +
  geom_bar(stat = 'identity') +
  ylim(0,100) +
  theme_minimal() +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/mC_tile_average_mC.pdf"), h = 2, w = 1.1)

write.table(plot_data_average, paste0(out_dir,"/", sample_to_crunch, "mC_tile_average_mC.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

###################
## Annotate domains

### CHG mC domains
mC_domains <- merged_mC_sans_orgs %>%
  mutate(domain_tmp = ifelse(chg_sites == "n", "no_sites", #this catches tiles with no CHG sites
                             ifelse(CHG >= MR_percent, "MR", #this catches sites with no data too: they get NA
                                    ifelse(CHG < UMR_percent, "UMR", "Intermediate")))) %>%
  mutate(CHG_based_domain = ifelse(is.na(domain_tmp), "Missing_Data", domain_tmp))

mC_domains

### 
# suggested filters
# UMR_percent 0.1
# MR_percent 0.4
# be careful this steps are maize specific and would need to be adjusted to another species...
# NOTE: include "& !is.na(CHH)" in any greater than tests or else will return NA rather than passing to next arg
# 2023/10/11 update - this starts with "merged_mC_sans_orgs" now and doesnt include the chg tile filter above
# 2023/10/11 update - fixed a few hundred occurances of CG-only tiles being called NA beccause of NA in chh or chg ratio columns...

mC_domains2 <- merged_mC_sans_orgs %>%
  mutate(domain_tmp = ifelse(chh_sites == "n", "no_sites", #this will catch the 1.45% of the genome that lack CHH sites (or any other site)
                      ifelse(CHH >= 0.15 & !is.na(CHH), "RdDM", #
                      ifelse(cg_sites == "n" & chg_sites == "y" & CHG >= MR_percent & !is.na(CHG), "Heterochromatin", # this catches tiles with CHG mC but no CG sites (otherwise they would get NA)
                      ifelse(cg_sites == "y" & chg_sites == "y" & CHG >= MR_percent & CG >= MR_percent  & !is.na(CG) & !is.na(CHG), "Heterochromatin",
                      ifelse(cg_sites == "y" & CG >= MR_percent & !is.na(CG), "CG_only",
                      ifelse(cg_sites == "y" & chg_sites == "y" & CG < UMR_percent & CHG < UMR_percent & CHH < UMR_percent & !is.na(CG) & !is.na(CHG) & !is.na(CHH), "Unmethylated",
                      ifelse(cg_sites == "n" & chg_sites == "y" & CHG < UMR_percent & CHH < UMR_percent & !is.na(CHG) & !is.na(CHH), "Unmethylated",
                      ifelse(cg_sites == "y" & chg_sites == "n" & CG < UMR_percent & CHH < UMR_percent & !is.na(CG) & !is.na(CHH), "Unmethylated",
                      ifelse(cg_sites == "y" & chg_sites == "y" & CG >= UMR_percent & CHG >= UMR_percent & CHH >= UMR_percent & !is.na(CG) & !is.na(CHG) & !is.na(CHH), "Intermediate",
                      ifelse(cg_sites == "n" & chg_sites == "y" & CHG >= UMR_percent & CHH >= UMR_percent & !is.na(CHG) & !is.na(CHH), "Intermediate",
                      ifelse(cg_sites == "y" & chg_sites == "n" & CG >= UMR_percent & CHH >= UMR_percent & !is.na(CG) & !is.na(CHH), "Intermediate",
                      ifelse(chg_sites == "y" & CHG >= UMR_percent & !is.na(CHG), "Heterochromatin", # this last one to catch some CHG-only tiles that were getting through and being called UMRs...
                      ifelse(cg_sites == "y" & CG >= UMR_percent & !is.na(CG), "Intermediate", # these are to catch any remaining intermediate mC, so its not counted as a UMR
                      ifelse(chg_sites == "y" & CHG >= UMR_percent & !is.na(CHG), "Intermediate", # these are to catch any remaining intermediate mC, so its not counted as a UMR
                      ifelse(chh_sites == "y" & CHH >= UMR_percent & !is.na(CHH), "Intermediate", # these are to catch any remaining intermediate mC, so its not counted as a UMR
                      ifelse(cg_sites == "n" | chg_sites == "n", "no_sites", NA))))))))))))))))) %>%
  mutate(domain = ifelse(is.na(domain_tmp), "Missing_Data", domain_tmp)) %>%
  mutate(domain_simple = ifelse(domain == "Heterochromatin", "MR",
                         ifelse(domain == "Missing_Data", "no_data",
                         ifelse(domain == "Unmethylated", "UMR",
                         ifelse(domain == "no_sites", "No_sites","other_mC")))))

mC_domains2

############
## summarise

mC_domains_freq <- mC_domains2 %>%
  group_by(domain) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  mutate(percent = total/sum(total)*100,
         MB = total*100/1000000)
mC_domains_freq

write.table(mC_domains_freq, paste0(out_dir,"/", sample_to_crunch, annotation_suffix,
                                    "_freq.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

# plot
domain_order <- c("Heterochromatin", "RdDM", "CG_only", "Unmethylated", "Intermediate", "no_sites", "Missing_Data")

mC_domains_freq_plot <- mC_domains_freq %>%
  gather(key = metric, value = number, -domain) %>%
  mutate(domain = factor(domain, levels = domain_order))

g <- ggplot(mC_domains_freq_plot, aes(x = domain, y = number)) +
  geom_bar(stat = 'identity') +
  facet_grid(metric ~ ., scales = 'free') +
  theme_minimal() +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/mC_domain_bar.pdf"), h = 4, w = 2)

######
mC_domains_freq_simple <- mC_domains2 %>%
  group_by(domain_simple) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  mutate(percent = total/sum(total)*100,
         MB = total*100/1000000)
mC_domains_freq_simple

write.table(mC_domains_freq_simple, paste0(out_dir,"/", sample_to_crunch, "_mC_domains_simple",
                                           "_cov_",coverage_filter_min,
                                           "_sites_",site_filter_min,
                                           "_MR_",MR_percent,
                                           "_UMR_",UMR_percent,
                                           "_freq.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

#mC_domains_freq_CHG <- mC_domains2 %>%
#  group_by(CHG_based_domain) %>%
#  summarise(total = n()) %>%
#  ungroup() %>%
#  mutate(percent = total/sum(total)*100,
#         MB = total*100/1000000)
#mC_domains_freq_CHG
#
#write.table(mC_domains_freq_CHG, paste0(out_dir,"/", sample_to_crunch, "_mC_domains_CHG",
#                                        "_cov_",coverage_filter_min,
#                                        "_sites_",site_filter_min,
#                                        "_MR_",MR_percent,
#                                        "_UMR_",UMR_percent,
#                                        "_freq.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)


######################### #########################

######################### #########################
# write bed file
mC_domains_bed <- mC_domains2 %>%
  mutate(start = start-1,
         score = ".",
         strand = ".") %>%
  select(chr, start, end, domain, score, strand)

write.table(mC_domains_bed, paste0(out_dir_beds,"/", sample_to_crunch, "_mC_domains",
                                   "_cov_",coverage_filter_min,
                                   "_sites_",site_filter_min,
                                   "_MR_",MR_percent,
                                   "_UMR_",UMR_percent,".bed"), sep = "\t", quote = F, row.names = F, col.names = F)

# write the whole data file
write.table(mC_domains2, paste0(out_dir_beds,"/", sample_to_crunch, "_mC_domains",
                                "_cov_",coverage_filter_min,
                                "_sites_",site_filter_min,
                                "_MR_",MR_percent,
                                "_UMR_",UMR_percent,".txt"), sep = "\t", quote = F, row.names = F, col.names = T)


## Make UMT and ND only bed files
# make UMR only bedfile
mC_domains2 %>% distinct(chr)

############# UMTs
# subset to UMTs
mC_domains2 %>% distinct(domain)
UMT_only <- mC_domains2 %>%
  filter(domain == "Unmethylated") %>%
  mutate(start = start-1) %>%
  select(chr:end, domain)
UMT_only
# 1,071,742

write.table(UMT_only, paste0(out_dir_beds, "/",sample_to_crunch, "_UMTs.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

############# Missing data and no sites
# subset to NDs
ND_only <- mC_domains2 %>%
  filter(domain %in% c("Missing_Data", "no_sites")) %>%
  mutate(start = start-1) %>%
  select(chr:end, domain)
ND_only
# 3,372,705

write.table(ND_only, paste0(out_dir_beds, "/", sample_to_crunch, "_NDs.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

############# Tiles with data (exclude no data and no sites)

NDs_data <- mC_domains_bed %>%
filter(!domain %in% c("no_sites", "Missing_Data")) %>%
select(chr, start, end, domain)

write.table(NDs_data, paste0(out_dir_beds,"/", sample_to_crunch, "_mC_domains",
                                   "_cov_",coverage_filter_min,
                                   "_sites_",site_filter_min,
                                   "_MR_",MR_percent,
                                   "_UMR_",UMR_percent, "_tiles_with_data",".bed"), sep = "\t", quote = F, row.names = F, col.names = F)
