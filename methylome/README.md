# Methylome pipeline

This is a multi purpose plant methylome pipeline that is used for processing SeqCap-v2 data, WGBS data, EM-seq data and also for identifying Unmethylated Regions (UMRs).

**Disclaimer:** these pipelines document the reproducible steps we take to generate the data on our server but may not be a plug and play on a different server. If you would like to reproduce the pipeline, you likely need to modify for your purposes or get in contact and we may be able to help.

--------

## Installation Notes
This pipeline requires a number of pieces of software, some we load as modules, some are expected to be in your path, depending on the step.

It is designed to run on UQ compute clusters.

Bsmap is a special case. It looks for bsmap in ~/software/bsmap-2.74 and puts this dir in your path. For historical reasons, the version of samtools shipped with Bsmap is used for most steps.


To do things:
 - soft code more steps (some things are hard coded - our bad)
 - described dependancies and required files structures

## Example pipeline execution - Server Steps

### Step 1 Trim reads

Use trimgalore to trim the reads.

usage="USAGE:
bash 01-trim_galore_qsub.sh <sample_list.txt>"

```
#01-trim_galore
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/01-trim_galore_qsub.sh \
samples.txt
```
**Methods summary**
Reads were trimmed and QC'ed with trim_galore version 0.4.3, powered by cutadapt v1.8.1 and fastqc v0.11.5.

### Step 2 Mapping using bsmap

Map using bsmap. This is preferred over bismark for speed with large genomes.

usage="USAGE:
bash 02-bsmap_qsub.sh <sample_list.txt> <genome.fa>"

```
#02-bsmap
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/02-bsmap_qsub.sh \
samples.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa
```

**Methods summary**
Reads were aligned with bsmap v2.74 with the following parameters -v 5 to allow allow 5 mismatches, -r 0 to report only unique mapping pairs, -p 1, -q 20 to allow quality trimming to q20, -A AGATCGGAAGAGC adapter sequence. Output file is in SAM format to allow custom QC and sorting (because some reads are not properly paired).

### Step 3 fix sam files

The output sam files produced by bsmap contain incorrect sam flags. Where PE reads map to different choromsomes the reads are still marked as correctly paird. This breaks picard tools, therefore, fix using samtools fixsam. Also make bams, sort and index (incase we want to vies these in IGV).

usage="USAGE:
bash 03-fix-sort-bsmap_qsub.sh <sample_list.txt>
"

```
#03-fix-sort-bsmap
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/03-fix-sort-bsmap_qsub.sh \
samples.txt
```

### Step 4 filter

Use picard tools to filter reads, removing duplicates, removing improperly paired reads and trimming overlapping reads so as to only count 'Cs' once per sequenced molecule. Also collect On target metrics.

usage="USAGE:
bash 04-filter_qsub.sh <sample_list.txt> <genome.fa> <CalculateHsMetrics_reference.bed>
for example:
"

```
#04-filter
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/04-filter_qsub.sh \
samples.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa \
/home/springer/pcrisp/ws/refseqs/maize/seqcapv2_onTarget-for-picard.bed
```

In the WGBS version the collect HS-metrics step is omited

```
#04-filter
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/04-filter-WGBS_qsub.sh \
samples.txt

```



### Step 5 Extract methylation data

Use methylratio.py script from bsmap to extract methylation data. Then use awk to convert output to mC context (CG, CHG, CHH) and parse coordinate to zero-based format "BED" type for bedtools. Then use awk to convert to bedgraph format (chr, start, stop, ratio) and split into a separate file per context. Use bedGraphToBigWig to make a bigWig file for IGV. Use awk to get conversion rates using the chloroplast reads. Use bedtools to intersect bam with target regions then use awk to sum reads per region. Then use bedtools to intersect C and CT counts file with target regions and use awk to sum counts. Also count methylaiton per region using eff_CT incase we want this metric later.

usage="USAGE:
bash 05-summarise_methylation_qsub.sh <sample_list.txt> <genome.fa> <intersect_regions.bed>
for example:
"

```
#05-summarise_methylation
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/05-summarise_methylation_qsub.sh \
samples.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa \
/home/springer/pcrisp/ws/refseqs/maize/BSseqcapv2_specific_regions.bed
```

output sites file: in folder BSMAPratio/...BSMAP_out.txt

chr     start end     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower        CI_upper
1       489     490     +       CHG     0.250   8.00    2       8       0       0       0.071   0.591
1       490     491     +       CG      0.875   8.00    7       8       0       0       0.529   0.978

output capture or target list summary: in folder OnTargetCoverage/...BSMAP_out_ontarget_mC.txt
Note: "sites" is sites with data, if C has no data it will not appear in the output, therefore not count towards the number of sites

<chr> <start> <stop> <context> <sites> <mC> <CT>
3       157881739       157881839       CHH     25      2       80
4       211836847       211836947       CHG     8       14      29

### Step 5.2 Context means

This step calculates means per context and sub-context. If a loci of interest file is provided it will additionally output a summary for these loci only. If there are no loci of interest then specify ``"none"``.

```
bash \
~/gitrepos/springerlab_methylation/SeqCap/05.2-analysis-contextMeans-qsub.sh \
samples.txt \
analysis/BSMAPratio \
20 \
none \
none
```

### Step 7 100 bp tiles

Summarise methylation to 100 bp tiles across the genome. This is a 4-parter. First it runs the perl script `met_context_window.pl` to get C, CT, ratios and C_site_counts per 100bp window per context. Per sample any window with data is reported for all contexts even if some windows do not have data for a particular context. Next the R script `07-tiles_bed_to_bigWig.R` is called to fix the ends of the chromosomes which have windows the extend beyond the chromosome ends, this is necessary if we want to make bigWigs etc. This script also pulls in the C_sites per window information from the file `maize_v4_100pb_tiles_zBased_sites.txt`. This file was created by me, tiles windows across the genome, then merges with the file ` AGPv4_cg_sites_012017.bed` from Jawon c/o Qing c/o Jackie which gives C sites per 100 bp window. Sites is actually a misleading term, this means individual C's per strand per context, is a CG site would count for 2 if both C's are in the window, a smaller fraction (1%) of CG sites bridge tiles, therefore meaning a small but significant number of tiles have 1,3,5 etc "C_sites". Note this C sites file does not include the contigs - might consider remaking this file to include the contigs (although they only account for 1.31% of the genome). Output is a bedgraph (fixed.sorted.bg; "chr", "start_zBased", "end", "ratio") and a text file with one-based coordinates (fixed.sorted.txt; "chr", "start", "end", "C", "CT", "ratio", "sites_with_data", "c_sites"). Then `bedGraphToBigWig` is used to make bigWigs from the tiles.

 usage="USAGE:
 bash 07-summarise_methylation_qsub.sh <sample_list.txt> <chrom.sizes file>
 for example:
 "

 Requires:
  - /home/springer/pcrisp/ws/refseqs/maize/maize_v4_100pb_tiles_zBased_sites.txt

```
#07-tiles_bed_to_bigWig
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/07-tiles_bed_to_bigWig_qsub.sh \
single_sample.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes
```

### 08-tiles_CHH_cov

This step takes the CHH output file from the 100 bp tiles script `*_BSMAP_out.txt.100.CHH.fixed.sorted.txt` and

1. Calculates CHH_cov (CHH_cov= CT/chh_sites).
2. Then using this field it filters on the arg $coverage_filter ($3); for example 2 or 5; which corrosponds to an average read coverage of 2 (or 5) accross all CHH sites in each 100 bp window, the output is significantly smaller file `*_BSMAP_out.txt.100.CHH_cov.txt`.

The idea (still being developed) is that:
1. These files can be used to explore the read coverage over the whole experiment; and
1. The tile index in this file could be combined with other samples to make a master list of tiles of interest (this is less memory intensive than combining the whole files).

usage="USAGE:
bash 08-tiles_analysis_qsub.sh <sample_list.txt> <data_folder> <coverage_filter>

For example:
```
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/08-tiles_analysis_CHH_cov_qsub.sh \
samples.txt \
tiles \
2
```
### 09-clean-up

Purges unnecesay files with the option to copy important files to home files space from scratch and also option to copy files to s3. Each of these three options (purge, copy to home and copy to s3) can be run independently.

See the script for more info and or the notebook.

### 10-tiles_filter_list

***Untested***

This step filters the `*_BSMAP_out.txt.100.CG.fixed.sorted.txt`, `*_BSMAP_out.txt.100.CHG.fixed.sorted.txt` and `*_BSMAP_out.txt.100.CHH.fixed.sorted.txt` files based on a list of tiles of interest. The idea is that once these smaller files are generated, they can be read into R in a stack for analysis all together.

For example:
```
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/10-tiles_analysis_filter_list_qsub.sh \
../samples.txt \
tiles \
../analysis_02_tiles_SeqCap_meta_140_samples/chh_2x_cov_62_sample_tile_list.tsv
```

--------

## Example pipeline execution - Laptop Steps

*It is possible these steps will be too intensive for the laptop with 100s of samples*

### Step 6 summarise-output

Use R to calculate ratios per region; add v4 annotations per region; and also merge with the read depth count file.

usage="USAGE:
bash 06-summarise-output-runner.sh <sample_list> <OnTargetCoverage_folder> <threads>
"

```
#06-summarise-output
bash ~/gitrepos/springerlab_methylation/SeqCap/06-summarise-output-runner.sh \
samples.txt \
OnTargetCoverage \
8
```

### Step 6b aggregate summarise-output

Use R to aggregate the summary output files into a single meta-table with all samples. Use this file to read into R and do further data exploration and calculate DMRs etc.

usage="USAGE:
Rscript 06-summarise-output-runner.sh <data_path> <outPrefix> <SeqCapEpi2_regions_annotation_v2_v4.csv_file>
"

```
Rscript \
~/gitrepos/springerlab_methylation/SeqCap/06b-aggregate-samples.R \
OnTargetCoverage_annotated \
Mei_final \
~/umn/refseqs/maize/SeqCap/Seqcap_ultimate_annotation_files/SeqCapEpi2_regions_annotation_v2_v4.csv
```

--------

# Scraping the log files

*Scrapers written so far*

1. Trim_galore - Parse the trimming logs to get the total reads.
2. Bsmap - to get initial mapping rates
2. MarkDuplicates - to get duplication rates
3. OnTargetMetrics - to get a heap of OnTargetMetrics but mainly, % On target, % Target 2X, mean target coverage, fold enrichment
4. ConversionRate - parse the output to get conversion rate
3. Methratio - get final valid mapping used to extract methylation data eg use to get total valid reads, C coverage, and use valid reads divided by total reads from the trim_galore logs to get mapping rate
4. CHH_cov - tiles passing 2x CHH

## Trimming (Total reads)

Parse the trimming logs to get the total reads. This is the total reads and should match PF clusters.

NOTE this will not match total reads passed to bsmap, because trimming will remove some reads entirely.

```
cd logs/..._01-trim_galore

for i in $(ls 01-trim_galore_e*); do
SAMPLE=$(grep '+ ID=' $i | cut -d "=" -f 2)
TOTAL_READS=$(grep 'Total number of sequences analysed:' $i | cut -d " " -f 6)
echo -e "$SAMPLE\t$TOTAL_READS"
done > total_reads_summary.tsv

d logs/..._01-trim_galore_swift

# updated with tr and squeeze to accommodate millions of reads samples
echo -e "Sample\tTotal_sequences_analysed\tPERCENT_READS_WITH_ADAPTERS_R1\tPERCENT_READS_WITH_ADAPTERS_R2\tPERCENT_BP_TRIMMED_R1\tPERCENT_BP_TRIMMED_R2" > ../total_reads_summary.tsv

for i in $(ls 01-trim_galore_swift_e*); do
SAMPLE=$(grep '+ ID=' $i | cut -d "=" -f 2)
TOTAL_READS=$(grep 'Total number of sequences analysed:' $i | tr -s ' ' | cut -d " " -f 6)
PERCENT_READS_WITH_ADAPTERS=$(grep 'Reads with adapters:' $i | tr -s ' ' | cut -d " " -f 5 | paste -sd '\t')
PERCENT_BP_TRIMMED=$(grep 'Quality-trimmed:' $i | tr -s ' ' | cut -d " " -f 4 | paste -sd '\t')
echo -e "$SAMPLE\t$TOTAL_READS\t$PERCENT_READS_WITH_ADAPTERS\t$PERCENT_BP_TRIMMED"
done >> ../total_reads_summary.tsv

cat ../total_reads_summary.tsv | column -t


# updated with tr and squeeze to accommodate millions of reads samples
# ATAC-seq .gz version
echo -e "Sample\tTotal_sequences_analysed\tPERCENT_READS_WITH_ADAPTERS_R1\tPERCENT_READS_WITH_ADAPTERS_R2\tPERCENT_BP_TRIMMED_R1\tPERCENT_BP_TRIMMED_R2" > ../total_reads_summary.tsv

for i in $(ls 01-trim_galore_gz_e*); do
SAMPLE=$(grep '+ ID=' $i | cut -d "=" -f 2)
TOTAL_READS=$(grep 'Total number of sequences analysed:' $i | tr -s ' ' | cut -d " " -f 6)
PERCENT_READS_WITH_ADAPTERS=$(grep 'Reads with adapters:' $i | tr -s ' ' | cut -d " " -f 5 | paste -sd '\t')
PERCENT_BP_TRIMMED=$(grep 'Quality-trimmed:' $i | tr -s ' ' | cut -d " " -f 4 | paste -sd '\t')
echo -e "$SAMPLE\t$TOTAL_READS\t$PERCENT_READS_WITH_ADAPTERS\t$PERCENT_BP_TRIMMED"
done >> ../total_reads_summary.tsv

cat ../total_reads_summary.tsv | column -t

```

## Bsmap Scraper

```
cd logs/..._02-bsmap

for i in $(ls 02-bsmap_e*); do
SAMPLE=$(grep 'echo sample being mapped is' $i | cut -d " " -f 7)
MAPPED_PAIRS=$(grep 'pairs:' $i | cut -d " " -f 8)
MAPPED_PAIRS_PERCENTAGE=$(grep 'pairs:' $i | cut -d " " -f 9 | cut -d "(" -f2 | cut -d ")" -f1)
MAPPED_A=$(grep 'single a:' $i | cut -d " " -f 6)
MAPPED_A_PERCENTAGE=$(grep 'single a:' $i | cut -d " " -f 7 | cut -d "(" -f2 | cut -d ")" -f1)
MAPPED_B=$(grep 'single b:' $i | cut -d " " -f 6)
MAPPED_B_PERCENTAGE=$(grep 'single b:' $i | cut -d " " -f 7 | cut -d "(" -f2 | cut -d ")" -f1)
echo -e "$SAMPLE\t$MAPPED_PAIRS\t$MAPPED_PAIRS_PERCENTAGE\t$MAPPED_A\t$MAPPED_A_PERCENTAGE\t$MAPPED_B\t$MAPPED_B_PERCENTAGE"
done > bsmap_summary.tsv

```

Total number of sequences analysed:

## BWA-meth

```
for i in $(ls 02-bwa-meth_o*); do \
SAMPLE=$(grep 'sample being mapped is' $i | tr -s ' ' | cut -d " " -f 5); \
MAPPED_READS=$(grep 'mapped (' $i); echo -e "$SAMPLE\t$MAPPED_READS"; \
done > bwa_summary.tsv

lst bwa_summary.tsv

```

## MarkDuplicates

MarkDuplicates - to get duplication rates

This step requires a samples.txt in the parent folder

```
# if before cleanup
cd /.../analysis/bsmapped_filtered
# if after cleanup
cd /.../analysis/HsMetrics_deDups_logs

for i in $(cat ../../samples.txt); do
SAMPLE=$i
UNPAIRED_READS_EXAMINED=$(grep 'Unknown Library' ${i}_MarkDupMetrics.txt | cut -f 2)
STATS=$(grep 'Unknown Library' ${i}_MarkDupMetrics.txt)
UNPAIRED_READS_EXAMINED=$(grep 'Unknown Library' ${i}_MarkDupMetrics.txt | cut -f 2)
echo -e "$SAMPLE\t$STATS"
done > MarkDuplicates_scraped.tsv
# add headder
echo -e 'sample\tLIBRARY\tBlurg\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE' \
| cat - MarkDuplicates_scraped.tsv > temp && mv temp MarkDuplicates_scraped.tsv

lst MarkDuplicates_scraped.tsv

```
## Overlapping Bases

```

echo -e "Sample\tPairs_overlapped\tAVERAGE_bp_Overlap" > ../overlapping_reads_summary.tsv

for i in $(ls 04-filter-WGBS-regular_e*); do
SAMPLE=$(grep '+ ID=' $i | cut -d "=" -f 2)
Pairs_overlapped=$(grep 'Number of overlapping pairs' $i | tr -s ' ' | cut -d " " -f 5)
AVERAGE_bp_Overlap=$(grep 'Reference Bases Overlapped' $i | tr -s ' ' | cut -d " " -f 6)
echo -e "$SAMPLE\t$Pairs_overlapped\t$AVERAGE_bp_Overlap"
done >> ../overlapping_reads_summary.tsv

cat ../overlapping_reads_summary.tsv | column -t

```

## OnTargetMetrics

```
cd ../analysis/bsmapped_filtered

#Get header from a random sample
HEADER=$(grep 'BAIT_SET' US_1_Index9_S14_HsMetrics_noDuplicate.txt | head -1)
echo -e "Sample\t$HEADER" > OnTargetMetrics_scraped.tsv

#grep 'BAIT_SET' *_HsMetrics_noDuplicate.txt | head -1 > OnTargetMetrics_scraped.tsv

#scrape
for i in $(cat ../../samples.txt); do
SAMPLE=$i
STATS=$(grep '^SeqCapEpi2_v4_capture_space_sorted' ${i}_HsMetrics_noDuplicate.txt)
echo -e "$SAMPLE\t$STATS"
done >> OnTargetMetrics_scraped.tsv

```

## ConversionRate

```
cd ../analysis/ConversionRate

echo -e "Sample\tC_counts\tCT_counts\tConversionRate" > ConversionRate_scraped.tsv

for i in $(cat ../../samples.txt); do
SAMPLE=$i
File_data=${i}_conversion_rate.txt
#STATS=$(cat ${i}_conversion_rate.txt)
C_counts=$(cat $File_data | cut -d " " -f 1)
CT_counts=$(cat $File_data | cut -d " " -f 2)
ConversionRate=$(cat $File_data | cut -d " " -f 3)
echo -e "$SAMPLE\t$C_counts\t$CT_counts\t$ConversionRate"
done >> ConversionRate_scraped.tsv
```


## Methratio

VALID_MAPPINGS, COVERED_CYTOSINES and AVERAGE_COVERAGE determined by methratio.py

Note sure how coverage is calculated...

```
cd logs/..._05-summarise_methylation

echo -e "sample\tVALID_MAPPINGS\tCOVERED_CYTOSINES\tAVERAGE_COVERAGE" > methratio_summary.tsv

for i in $(ls 05-summarise_methylation_o*); do
SAMPLE=$(grep 'sample being mapped is' $i | cut -d " " -f 5)
VALID_MAPPINGS=$(grep 'total' $i | cut -d " " -f 2)
COVERED_CYTOSINES=$(grep 'total' $i | cut -d " " -f 5)
AVERAGE_COVERAGE=$(grep 'total' $i | cut -d " " -f 10)
echo -e "$SAMPLE\t$VALID_MAPPINGS\t$COVERED_CYTOSINES\t$AVERAGE_COVERAGE"
done >> methratio_summary.tsv

# or WGBS

cd logs/..._05-summarise_methylation-WGBS

echo -e "sample\tVALID_MAPPINGS\tCOVERED_CYTOSINES\tAVERAGE_COVERAGE" > methratio_summary.tsv

for i in $(ls 05-summarise_methylation-WGBS_o*); do
SAMPLE=$(grep 'sample being mapped is' $i | cut -d " " -f 5)
VALID_MAPPINGS=$(grep 'total' $i | cut -d " " -f 2)
COVERED_CYTOSINES=$(grep 'total' $i | cut -d " " -f 5)
AVERAGE_COVERAGE=$(grep 'total' $i | cut -d " " -f 10)
echo -e "$SAMPLE\t$VALID_MAPPINGS\t$COVERED_CYTOSINES\t$AVERAGE_COVERAGE"
done >> methratio_summary.tsv
```

## CHH_cov

```
cd logs/..._08-tiles_analysis_CHH_cov

echo -e "sample\tTILES_CHH_COV" > CHH_cov_summary.tsv

for i in $(ls 08-tiles_analysis_CHH_cov_o*); do
SAMPLE=$(grep 'sample being mapped is' $i | cut -d " " -f 5)
TILES_CHH_COV=$(grep '\[1\] "For sample*' $i | cut -d " " -f 7)
echo -e "$SAMPLE\t$TILES_CHH_COV"
done >> CHH_cov_summary.tsv

```
---

# Section 2 - calling UMRs

This section describes the steps and input/output files for UMR calling.

**Disclaimer:** as above - this pipeline documents the reproducible steps we take to generate the data on our server including submission scripts etc but is unlikely to be a plug and play on a different server. If you would like to reproduce the pipeline, you likely need to modify for your purposes or get in contact and we may be able to help.

**Steps:**
1. Create a reference file with the coordinates of each 100 bp tile in the genome

2. Summarise per cytosine methylation calls (eg the output of bsmap-2.74/methratio.py) into average methylation per CG, CHG and CHH per 100 bp tile of the genome.

2. Categorise each tile into methylation domains (six domains or types, including “missing data” (“no data” and “no sites”), “High CHH/RdDM”, “Heterochromatin”, “CG only”, “Unmethylated” or “intermediate”, in preferential order as per Crisp et al. (2020) https://doi.org/10.1073/pnas.2010250117) and merge adjacent tiles with the same methylation type, then keep the UMR tiles.

## Step 1 - 100bp tile reference files

### code step - create_genome_tile_mC_counts

This step generates a file with the coordinates of 100bp tiles in the genome and the count of CG, CHG and CHH sites in each tile.

Dependencies (what we use):
- R v3.5.0
- bedtools v2.25.0
- seqkit v0.14.0
- csvtk v0.21.0
- the associated `create_genome_tile_mC_counts` R scripts and the packages specified in the R scripts
- seqtk 1.3-r115-dirty

---
Input
- a genome fasta file with the suffix `'fa'` eg `Sbicolor_454_v3.0.1.fa`

- a text file with the prefix of the genome to be indexed eg:
```
cat genome.txt
> Sbicolor_454_v3.0.1
>
(* on some systems a trailing empty line is needed if submitting a single job as a batch job)
```

- a text file with the chromosome sizes, with the same prefix and the suffix chrom.sizes - this can be generated with seqtk comp; eg:
```
seqtk comp Sbicolor_454_v3.0.1.fa | cut -f 1-2 > Sbicolor_454_v3.0.1.chrom.sizes
# <chromosome> <size>
cat Sbicolor_454_v3.0.1.chrom.sizes
> Chr01   80884392
> Chr02   77742459
> Chr03   74386277
> Chr04   68658214
> Chr05   71854669
> Chr06   61277060
> Chr07   65505356
> Chr08   62686529
> Chr09   59416394
> Chr10   61233695
```
---
Then run the `create_genome_tile_mC_counts-qsub.sh` script in the directory that contains all of these files

```
#create_genome_tile_mC_counts
bash \
/home/uqpcrisp/gitrepos/crisplab_epigenomics/methylome/create_genome_tile_mC_counts-qsub.sh \
genome.txt \
12:00:00 \
40
```
---

The key output is file named `${genome}_100bp_tiles_zBased_sites_counts.txt` this is the file you need for the next step and has the tiles and site counts.

```
# Sbicolor_454_v3.0.1_100bp_tiles_zBased_sites_counts.txt
chr     start   end     start_zBased    cg_sites        chg_sites       chh_sites
Chr01   1       100     0       2       2       42
Chr01   101     200     100     2       2       29
Chr01   201     300     200     2       2       12
Chr01   301     400     300     1       0       28
Chr01   401     500     400     5       0       26
Chr01   501     600     500     0       0       34
Chr01   601     700     600     6       6       19
Chr01   701     800     700     2       5       27
Chr01   801     900     800     2       9       32
Chr01   901     1000    900     8       7       35
Chr01   1001    1100    1000    6       2       26
```


## Step 2 Summarise methylation into 100 bp tiles

### code step - 07-tiles_bed_to_bigWig

Summarise methylation to 100 bp tiles across the genome. This is a 4-parter.

Dependencies (what we use):
- perl v5.16.3
- R v3.5.0
- 07-tiles_bed_to_bigWig.R
- bedGraphToBigWig (UCSC)

---

The input files are are generated by the previous step, or these could be generated by other means, if in the same format
- A text file with sample name prefix(s)
- the chromomsome sizes files used in the last step, in our example: `Sbicolor_454_v3.0.1.chrom.sizes`
- the `${genome}_100bp_tiles_zBased_sites_counts.txt` file generated in the last step, in our example: `Sbicolor_454_v3.0.1_100bp_tiles_zBased_sites_counts.txt`
- a file called `${genome}_BSMAP_out.txt`, this is the parsed output of `bsmap-2.74/methratio.py` in zero based format no header ie bed-like; example from wheat with a header (remove header):
```
chr          start  end  strand  context  ratio  eff_CT_count  C_count  CT_count  rev_G_count  rev_GA_count  CI_lower  CI_upper
chr4B_part1  25     26   +       CHH      0.000  21.00         0        21        0            0             0.000     0.155
chr4B_part1  31     32   +       CHH      0.087  23.00         2        23        0            0             0.024     0.268
chr4B_part1  45     46   +       CHH      0.000  30.00         0        30        0            0             -0.000    0.114
chr4B_part1  56     57   +       CG       0.097  31.00         3        31        0            0             0.033     0.249
chr4B_part1  66     67   +       CHH      0.067  30.00         2        30        0            0             0.018     0.213
chr4B_part1  78     79   +       CHH      0.031  32.00         1        32        0            0             0.006     0.157
chr4B_part1  80     81   +       CHH      0.000  36.00         0        36        0            0             0.000     0.096
chr4B_part1  84     85   +       CHG      0.029  35.00         1        35        0            0             0.005     0.145
chr4B_part1  85     86   +       CG       0.108  37.00         4        37        0            0             0.043     0.247
chr4B_part1  90     91   +       CHH      0.000  39.00         0        39        0            0             0.000     0.090
```
---

First this step runs the perl script `met_context_window.pl` to get C, CT, ratios a summarised per 100bp window per context. Per sample any window with data is reported for all contexts even if some windows do not have data for a particular context. Next the R script `07-tiles_bed_to_bigWig.R` is called to fix the ends of the chromosomes which have windows that extend beyond the chromosome ends, this is necessary if we want to make bigWigs etc. This script also pulls in the C_sites per window information from the file `${genome}_100bp_tiles_zBased_sites_counts.txt` to enable filter out tiles with less than the specified threshold of "Cs" per context, as desired.

```
 usage="USAGE:
 bash 07-summarise_methylation_qsub.sh <sample_list.txt> <chrom.sizes file> <reference tile file> <walltime> <cores>
 for example:
 "
```
In our work flow this is run from the working directory of the project and expects to find a subdirectory layout like below, the output is written to the folder called `tiles`

```
example_project/
├── analysis
│   ├── BSMAPratio
│   │   └── Chinese_Spring_BSMAP_out.txt
│   └── tiles
│       ├── Chinese_Spring_BSMAP_out.txt.100.CG.fixed.sorted.txt
│       ├── Chinese_Spring_BSMAP_out.txt.100.CHG.fixed.sorted.txt
│       └── Chinese_Spring_BSMAP_out.txt.100.CHH.fixed.sorted.txt
├── logs
└── samples.txt
```

```
cat samples.txt
> Chinese_Spring
> Awesome_variety_2
> etc
```

Example command to call the script

```
#07-tiles_bed_to_bigWig
bash \
/home/uqpcrisp/gitrepos/crisplab_epigenomics/methylome/07-tiles_bed_to_bigWig_qsub.sh \
samples.txt \
/90days/uqpcrisp/tmp_refseqs/sorghum/Sbicolor_454_v3.0.1.chrom.sizes \
/90days/uqpcrisp/tmp_refseqs/sorghum/Sbicolor_454_v3.0.1_100bp_tiles_zBased_sites_counts.txt \
4:00:00 \
30
```


Output is ... (for each sample processed)

```
# Chinese_Spring_BSMAP_out.txt.100.CG.fixed.sorted.txt
chr          start      end        C     CT    ratio                sites_with_data  cg_sites
chr1A_part1  101        200        38    40    0.95                 12               12
chr1A_part1  201        300        88    108   0.814814814814815    16               16
chr1A_part1  301        400        235   290   0.810344827586207    30               30
chr1A_part1  401        500        127   136   0.933823529411765    12               12
chr1A_part1  501        600        281   321   0.875389408099688    14               14
chr1A_part1  601        700        101   146   0.691780821917808    10               10
chr1A_part1  701        800        26    30    0.866666666666667    8                8
chr1A_part1  801        900        15    19    0.789473684210526    13               14
```

```
# Chinese_Spring_BSMAP_out.txt.100.CHG.fixed.sorted.txt
chr          start      end        C     CT    ratio                sites_with_data  chg_sites
chr1A_part1  101        200        15    25    0.6                  9                10
chr1A_part1  201        300        44    77    0.571428571428571    11               11
chr1A_part1  301        400        71    150   0.473333333333333    14               14
chr1A_part1  401        500        116   185   0.627027027027027    13               13
chr1A_part1  501        600        2     16    0.125                2                2
chr1A_part1  601        700        101   206   0.490291262135922    14               14
chr1A_part1  701        800        5     15    0.333333333333333    3                3
chr1A_part1  801        900        4     12    0.333333333333333    8                12
```

```
# Chinese_Spring_BSMAP_out.txt.100.CHH.fixed.sorted.txt
chr          start      end        C     CT    ratio                sites_with_data  chh_sites
chr1A_part1  1          100        0     2     0                    2                43
chr1A_part1  101        200        2     50    0.04                 21               29
chr1A_part1  201        300        20    217   0.0921658986175115   31               31
chr1A_part1  301        400        42    237   0.177215189873418    25               25
chr1A_part1  401        500        23    308   0.0746753246753247   22               22
chr1A_part1  501        600        43    646   0.0665634674922601   30               30
chr1A_part1  601        700        35    538   0.0650557620817844   29               29
chr1A_part1  701        800        6     140   0.0428571428571429   36               36
chr1A_part1  801        900        3     56    0.0535714285714286   26               31
chr1A_part1  1001       1100       0     5     0                    5                22
chr1A_part1  1101       1200       0     30    0                    28               30
```


# Step 3 UMR calling

The following 6 files are usually produced in UMR identification.

The following 3 files were processed requiring a minimum coverage of at least 3x per cytosine per strand.

Chinese_Spring_cov_3_sites_2_MR_0.4_UMR_0.1_UMRs_6col.bed

Chinese_Spring_mC_domains_cov_3_sites_2_MR_0.4_UMR_0.1_tiles_with_data.bed

Chinese_Spring_mC_domains_cov_3_sites_2_MR_0.4_UMR_0.1_tiles_with_data_inc_NDs_merged.bed



The following 3 files were processed requiring a minimum coverage of at least 5x per cytosine per strand.

Chinese_Spring_cov_5_sites_2_MR_0.4_UMR_0.1_UMRs_6col.bed

Chinese_Spring_mC_domains_cov_5_sites_2_MR_0.4_UMR_0.1_tiles_with_data.bed

Chinese_Spring_mC_domains_cov_5_sites_2_MR_0.4_UMR_0.1_tiles_with_data_inc_NDs_merged.bed

File Descriptions:

"*UMRs_6col.bed"
- this is the coordinates of the UMRs
- columns are: <chr> <start> <stop> <methylation_type> <UMR_ID> <size_category>
- methylation_type is always "UMR"
- size_category small (< 300), medium (< 900) or large (>=900)
- However, regions less than 300 bp are excluded from the final UMR classification.


"*_tiles_with_data.bed"
- this is the co-ordinates of all tiles with data at the specific threshold (3x or 5x in this example) and with the minimum number of cytosines (2 in this case)
- columns are: <chr> <start> <stop> <methylation_type>
- the definition of methylation type is described in the below heuristic:

The methylation categories include “missing data” (including “no data” and “no sites”), “RdDM,” “heterochromatin,” “CG-only,” “unmethylated,” or “intermediate”. In this analysis we are primarily interested in using this classification to identify the UMRs; however, users might also be interested in other types of methylation, which could be extracted from this same analysis strategy. We also suggest removing organelles from the data before proceeding with this step; however, in this example data the organelle genomes have already been removed.

This analysis is performed using a custom R script that we have provided Call-umrs.R. Regions are classified according in the following hierarchy: tiles are classified as missing data if tiles have less than two cytosines in the relevant context or if there is less than the specified coverage threshold of reads (eg 3-5× coverage); RdDM if CHH methylation is greater than 15%; heterochromatin if CG and CHG methylation is 40% or greater; CG-only if CG methylation is greater than 40%; unmethylated if CG, CHG, and CHH are less than 10%; and intermediate if methylation is 10% or greater but less than 40%. Note that the levels of CHH methylation are hard coded in this script; while the level of CG and CHG are specified when calling the script. We have found these levels to be appropriate for a range of species; however, they could be adjusted if your genome has a different or unusual distribution, for example if CHH methylation is known to be higher.


"_tiles_with_data_inc_NDs_merged.bed"
- this is the co-ordinates of all tiles used in the UMR classification including some tiles originally classified as "no data"
Following UMR tile classification, adjacent unmethylated tiles (UMTs) were merged. To capture and combine any unmethylated regions that were fragmented by a short interval of missing data (low coverage or no sites), any merged UMT regions that were separated by missing data were also merged so long as the resulting merged region consisted of no more than 33% missing data.


Dependencies (what we use):
- R v3.5.0
- bedtools v2.25.0
- the associated `21-call-umrs` R scripts and the packages specified in the R scripts

Input
- the chromomsome sizes files used in the last step, in our example: `Sbicolor_454_v3.0.1.chrom.sizes`
- the `${genome}_100bp_tiles_zBased_sites_counts.txt` file generated in the last step, in our example: `Sbicolor_454_v3.0.1_100bp_tiles_zBased_sites_counts.txt`

```
bash \
/home/uqpcrisp/gitrepos/crisplab_epigenomics/methylome/21-call-umrs-qsub.sh \
../samples.txt \
/90days/uqpcrisp/tmp_refseqs/sorghum/Sbicolor_454_v3.0.1_100bp_tiles_zBased_sites_counts.txt \
/90days/uqpcrisp/tmp_refseqs/sorghum/Sbicolor_454_v3.0.1_sorted.chrom.sizes \
5 \
2 \
0.4 \
0.1 \
00:30:00 \
20
```
