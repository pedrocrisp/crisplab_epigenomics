#!/bin/bash -l
#PBS -A UQ-SCI-SAFS
#PBS -N sum_meth-WGBS-wheat2
#PBS -r y
#PBS -m abej
#PBS -M p.crisp@uq.edu.au

# for barley increase walltime to 48hr + 80 Gb

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: array_ID is ${PBS_ARRAY_INDEX}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

########## Modules #################
# python 2.7.8 no longer available
# module load python2/2.7.8
module load python/2.7.14
#module load java
module load bedtools
#module load bamtools
#bsmap requires samtools < 1.0.0
# module load samtools/0.1.18
PATH=~/software/bsmap-2.74/samtools:$PATH

LC_COLLATE=C

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p BSMAPratio
mkdir -p TempOut
#mkdir -p OnTargetCoverage
mkdir -p ConversionRate

########## Run #################
        # The required input is all in the folder bsmapped_filtered from the prior step
        ########################
        # extract methylation information using bsmap tool methratio.py

#         python /home/uqpcrisp/software/bsmap-2.74/methratio.py \
#         -o BSMAPratio/${ID}_methratio_chr1.txt \
#         -d ${genome_reference} \
#         -u \
#         -z \
#         -c chr1A_part1,chr1A_part2,chr1B_part1,chr1B_part2,chr1D_part1,chr1D_part2 \
#         -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
#
# echo ch1 summarised
#
#         python /home/uqpcrisp/software/bsmap-2.74/methratio.py \
#         -o BSMAPratio/${ID}_methratio_chr2.txt \
#         -d ${genome_reference} \
#         -u \
#         -z \
#         -c chr2A_part1,chr2A_part2,chr2B_part1,chr2B_part2,chr2D_part1,chr2D_part2 \
#         -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
#
# echo ch2 summarised
#
#         python /home/uqpcrisp/software/bsmap-2.74/methratio.py \
#         -o BSMAPratio/${ID}_methratio_chr3.txt \
#         -d ${genome_reference} \
#         -u \
#         -z \
#         -c chr3A_part1,chr3A_part2,chr3B_part1,chr3B_part2,chr3D_part1,chr3D_part2 \
#         -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
#
# echo ch3 summarised
#
        python /home/uqpcrisp/software/bsmap-2.74/methratio.py \
        -o BSMAPratio/${ID}_methratio_chr4.txt \
        -d ${genome_reference} \
        -u \
        -z \
        -c chr4A_part1,chr4A_part2,chr4B_part1,chr4B_part2,chr4D_part1,chr4D_part2 \
        -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

echo ch4 summarised
#
#         python /home/uqpcrisp/software/bsmap-2.74/methratio.py \
#         -o BSMAPratio/${ID}_methratio_chr5.txt \
#         -d ${genome_reference} \
#         -u \
#         -z \
#         -c chr5A_part1,chr5A_part2,chr5B_part1,chr5B_part2,chr5D_part1,chr5D_part2 \
#         -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
#
# echo ch5 summarised
#
#         python /home/uqpcrisp/software/bsmap-2.74/methratio.py \
#         -o BSMAPratio/${ID}_methratio_chr6.txt \
#         -d ${genome_reference} \
#         -u \
#         -z \
#         -c chr6A_part1,chr6A_part2,chr6B_part1,chr6B_part2,chr6D_part1,chr6D_part2 \
#         -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
#
# echo ch6 summarised
#
#         python /home/uqpcrisp/software/bsmap-2.74/methratio.py \
#         -o BSMAPratio/${ID}_methratio_chr7.txt \
#         -d ${genome_reference} \
#         -u \
#         -z \
#         -c chr7A_part1,chr7A_part2,chr7B_part1,chr7B_part2,chr7D_part1,chr7D_part2,chrUn \
#         -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
#
# echo ch7 summarised
#
#         # grab header
#         head -n 1 BSMAPratio/${ID}_methratio_chr1.txt > BSMAPratio/${ID}_methratio.txt
#
#         # cat files without header
#         find BSMAPratio/ -type f -name "${ID}_methratio*" -print | while read filename; do
#             # echo ${filename}
#             tail -n +2 ${filename}
#         done >> BSMAPratio/${ID}_methratio.txt

# rename to comply with rest of pipeline -
# THIS ONLY WORKS IF YOU ONLY MADE A SINGLE SUMMARY FILE PER SAMPLE
# IF YOU RUN MULTIPLE CHR SEPERATELY, THEN CAT TOGTETHER FIRST
        mv BSMAPratio/${ID}_methratio* BSMAPratio/${ID}_methratio.txt

        #awk funciton for extracting methylation info from methratio.py output. Check with Qing what this is meant to do. Also try to figure out how to split this over multiple lines
        #awk '(NR>1){if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHH""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else print $1"\t"$2-1"\t"$2"\t"$3"\t""CNN""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' BSMAPratio/${ID}.txt > BSMAPratio/${ID}_BSMAP_out.txt

        # Well 2nt resolution makes a proper zero-based coordinate BED file... (this is what Qing did).
        awk_make_bed='BEGIN {OFS = FS} (NR>1){
                if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./))
                        print $1, $2-1, $2, $3, "CG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/))
                        print $1, $2-1, $2, $3, "CHG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/))
                        print $1, $2-1, $2, $3, "CHH", $5, $6, $7, $8, $9, $10, $11, $12;
                else
                        print $1, $2-1, $2, $3, "CNN", $5, $6, $7, $8, $9, $10, $11, $12
                }
                '

        # make subcontext bed
        # CG
        # CHG sub context and minus strand sequence (reverse complement):
        # CAG, CCG, CTG
        # CTG, CGG, GAC
        # CHH subcontext and reverse complement
        # CAA, CAC, CAT, CCA, CCC, CCT, CTA, CTC, CTT
        # TTG, GTG. ATG, TGG, GGG, AGG, TAG, GAG, AAG

        awk_make_subcontext_bed='BEGIN {OFS = FS} (NR>1){
                if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./))
                        print $1, $2-1, $2, $3, "CG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^CTG../ ) || ($3=="+" &&  $4~/^..CAG/))
                        print $1, $2-1, $2, $3, "CAG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^CGG../ ) || ($3=="+" &&  $4~/^..CCG/))
                        print $1, $2-1, $2, $3, "CCG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GAC../ ) || ($3=="+" &&  $4~/^..CTG/))
                        print $1, $2-1, $2, $3, "CTG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^TTG../ ) || ($3=="+" &&  $4~/^..CAA/))
                        print $1, $2-1, $2, $3, "CAA", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GTG../ ) || ($3=="+" &&  $4~/^..CAC/))
                        print $1, $2-1, $2, $3, "CAC", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^ATG../ ) || ($3=="+" &&  $4~/^..CAT/))
                        print $1, $2-1, $2, $3, "CAT", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^TGG../ ) || ($3=="+" &&  $4~/^..CCA/))
                        print $1, $2-1, $2, $3, "CCA", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GGG../ ) || ($3=="+" &&  $4~/^..CCC/))
                        print $1, $2-1, $2, $3, "CCC", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^AGG../ ) || ($3=="+" &&  $4~/^..CCT/))
                        print $1, $2-1, $2, $3, "CCT", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^TAG../ ) || ($3=="+" &&  $4~/^..CTA/))
                        print $1, $2-1, $2, $3, "CTA", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GAG../ ) || ($3=="+" &&  $4~/^..CTC/))
                        print $1, $2-1, $2, $3, "CTC", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^AAG../ ) || ($3=="+" &&  $4~/^..CTT/))
                        print $1, $2-1, $2, $3, "CTT", $5, $6, $7, $8, $9, $10, $11, $12;
                else
                        print $1, $2-1, $2, $3, "CNN", $5, $6, $7, $8, $9, $10, $11, $12
                }
                '

        #awk -F$"\t" "$awk_make_bed" "F1-16_Index5_S1_methratio.txt" > F1-16_Index5_S1.bed
        #BSMAPratio/$F1-16_Index5_S1_BSMAP_out.txt

        awk -F$"\\t" "$awk_make_bed" \
        "BSMAPratio/${ID}_methratio.txt" > "BSMAPratio/${ID}_BSMAP_out.txt"

        if [ "$make_subcontext" == "yes" ]
        then
        awk -F$"\\t" "$awk_make_subcontext_bed" \
        "BSMAPratio/${ID}_methratio.txt" > "BSMAPratio/${ID}_BSMAP_out_subcontext.txt"
        fi
        ########################
        #For genome browser

        # begGraph ratio files for tdfs
        awk_make_bedGraph='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $8/$9*100, $5
        }
        '

        # split bedGraph by contex
        awk_make_bedGraph_context='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $4 > "BSMAPratio/"ID"_BSMAP_out_"$5".bedGraph"
        }
        '

        if [ "$make_subcontext" == "yes" ]
        then
        # split bedGraph by sub-contex
        # only difference in the output file suffix (so that we make CG files for context and sub-contex: sainty check - they should be the same)
        awk_make_bedGraph_subcontext='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $4 > "BSMAPratio/"ID"_BSMAP_out_subcontext_"$5".bedGraph"
        }
        '
        fi

        #pipe bedGraph to split by context (use dash to read from sdtin)
        # per context
        awk -F$"\\t" "$awk_make_bedGraph" \
        "BSMAPratio/${ID}_BSMAP_out.txt" | \
        awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" -

        if [ "$make_subcontext" == "yes" ]
        then
        # per sub-context
        awk -F$"\\t" "$awk_make_bedGraph" \
        "BSMAPratio/${ID}_BSMAP_out_subcontext.txt" | \
        awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_subcontext" -
        fi

        # sort
        sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_CG.bedGraph > BSMAPratio/${ID}_BSMAP_out_CG.sorted.bedGraph
        sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_CHG.bedGraph > BSMAPratio/${ID}_BSMAP_out_CHG.sorted.bedGraph
        sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_CHH.bedGraph > BSMAPratio/${ID}_BSMAP_out_CHH.sorted.bedGraph

        #Make bigWigs per context
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CG.sorted.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_CG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CHG.sorted.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_CHG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CHH.sorted.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_CHH.bigWig"

        if [ "$make_subcontext" == "yes" ]
        then
        #Make bigWigs per CHG sub context
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CG.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CG.bigWig"

        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CAG.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CAG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CCG.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CCG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CTG.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CTG.bigWig"

        #Make bigWigs per CHH sub context
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CAA.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CAA.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CAC.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CAC.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CAT.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CAT.bigWig"

        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CCA.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CCA.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CCC.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CCC.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CCT.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CCT.bigWig"

        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CTA.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CTA.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CTC.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CTC.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_subcontext_CTT.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_subcontext_CTT.bigWig"
        fi

        #remove bedGraph it is large and not really required
        # keep bigWigs

        # keeing temporarily in case useful for metaplots? Does it have zeros?
        rm -rv BSMAPratio/${ID}*.bedGraph

        ########################
                ## Now make stranded bigWigs for IGV
                # make both sets, decided later which is prefered

                ## make begGraph ratio files for bigWigs
                # function
                awk_make_bedGraph='BEGIN {OFS = FS} (NR>1){
                  print $1, $2, $3, $8/$9*100, $5, $4
                }
                '
                # run function
                awk -F$"\\t" "$awk_make_bedGraph" \
                "BSMAPratio/${ID}_BSMAP_out.txt" > "BSMAPratio/${ID}_BSMAP_out.bg"

                ## split into + and - strand based on column 6
                awk -F$"\\t" -v ID=$ID 'BEGIN {OFS = FS} (NR>1){
                  print > "BSMAPratio/"ID"_BSMAP_out_"$6".bedGraph"
                }' "BSMAPratio/${ID}_BSMAP_out.bg"

                ## split plus and minnus bedGraphs by contex
                # plus funciton
                awk_make_bedGraph_context_plus='BEGIN {OFS = FS} (NR>1){
                  print $1, $2, $3, $4 > "BSMAPratio/"ID"_BSMAP_out_+_"$5".bedGraph"
                }
                '
                # minus funciton (multiply ratio by -1)
                awk_make_bedGraph_context_minus='BEGIN {OFS = FS} (NR>1){
                  print $1, $2, $3, $4*-1 > "BSMAPratio/"ID"_BSMAP_out_-_"$5".bedGraph"
                }
                '
                # run functions
                awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context_plus" "BSMAPratio/"${ID}"_BSMAP_out_+.bedGraph"
                awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context_minus" "BSMAPratio/"${ID}"_BSMAP_out_-.bedGraph"

                # sort
                sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_+_CG.bedGraph > BSMAPratio/${ID}_BSMAP_out_+_CG.sorted.bedGraph
                sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_+_CHG.bedGraph > BSMAPratio/${ID}_BSMAP_out_+_CHG.sorted.bedGraph
                sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_+_CHH.bedGraph > BSMAPratio/${ID}_BSMAP_out_+_CHH.sorted.bedGraph

                sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_-_CG.bedGraph > BSMAPratio/${ID}_BSMAP_out_-_CG.sorted.bedGraph
                sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_-_CHG.bedGraph > BSMAPratio/${ID}_BSMAP_out_-_CHG.sorted.bedGraph
                sort -k1,1 -k2,2n BSMAPratio/${ID}_BSMAP_out_-_CHH.bedGraph > BSMAPratio/${ID}_BSMAP_out_-_CHH.sorted.bedGraph


                ## Make bigWigs per context
                # bigWigs for each individual induction with its own genome reference and plasmid
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_+_CG.sorted.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_+_CG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_+_CHG.sorted.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_+_CHG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_+_CHH.sorted.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_+_CHH.bigWig"

                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_-_CG.sorted.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_-_CG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_-_CHG.sorted.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_-_CHG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_-_CHH.sorted.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_-_CHH.bigWig"

        #remove bedGraph it is large and not really required
        # keep bigWigs
        rm -rv BSMAPratio/${ID}*.bedGraph

        ########################

        # conversion rate
        #awk -F"\t" '{if($1=="Pt") print}' "./BSMAPratio/"${ID}"_BSMAP_out.txt" | awk '{sum1 += $8; sum2 +=$9} END {print sum1"\t"sum2"\t"100-sum1/sum2*100}' > "./ConversionRate/"$i"_conversion_rate.txt"

        # check if there is PT CHH data
        PT_data=$(awk -F$"\\t" \
        'BEGIN {OFS = FS} {if($1=="Pt" && $5=="CHH") print}' \
        BSMAPratio/${ID}_BSMAP_out.txt |wc -l)
        echo $PT_data

        # conversion rate pete - not sure if this is correct...
        # using plastid CHH unconverted C rate: 100-(sum(#C_counts)/sum(#CT_counts)*100)
        # Note that bsmap recomends using the eff_CT_counts (field #$7) which considers if there is a mismatch with the reverse strand.
        # However, Qing recommends just using the CT count (field #$9) becasue the reverse strand could equally be a sequencing error. Check this.
        if [ $PT_data -gt 0 ]
        then
          echo 'Calculating conversion rate'
          awk -F$"\\t" \
          'BEGIN {OFS = FS} {if($1=="Pt" && $5=="CHH") print}' \
          BSMAPratio/${ID}_BSMAP_out.txt | \
          awk '{sum1 += $9; sum2 +=$8} END {print sum1, sum2 , 100-((sum2/sum1)*100)}' > ConversionRate/${ID}_conversion_rate.txt
          # conversion rate pete - using eff_CT
          awk -F$"\\t" \
          'BEGIN {OFS = FS} {if($1=="Pt" && $5=="CHH") print}' \
          BSMAPratio/${ID}_BSMAP_out.txt | \
          awk '{sum1 += $7; sum2 +=$8} END {print sum1, sum2 , 100-((sum2/sum1)*100)}' > ConversionRate/${ID}_conversion_rate_eff_C.txt
        else
          echo 'No data for conversion rate calculation'
        fi
        #debugging
        #awk -F$"\\t" \
        #'BEGIN {OFS = FS} {if($1=="Pt" && $5=="CHH") print}' \
        #BSMAPratio/F1-16_Index5_S1_BSMAP_out.txt | \
        #awk -F$"\\t" 'BEGIN {OFS = FS} {sum1 += $7; sum2 +=$8} END {print sum1, sum2 , 100-((sum2/sum1)*100)}'


echo finished summarising
