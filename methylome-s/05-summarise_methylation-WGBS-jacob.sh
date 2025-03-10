#!/bin/bash -l
#SBATCH --job-name summarise_methylation-WGBS
#SBATCH --requeue
#SBATCH --partition=general

# for barley increase walltime to 48hr + 80 Gb

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo SBATCH: working directory is $SLURM_SUBMIT_DIR
echo SBATCH: job identifier is $SLURM_JOBID
echo SBATCH: array_ID is ${SLURM_ARRAY_TASK_ID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to SLURM_SUBMIT_DIR
cd "$SLURM_SUBMIT_DIR"
echo working dir is now $PWD

########## Modules #################
# python 2.7.8 no longer available
# module load python2/2.7.8
# module load python/2.7.14
module load python/2.7.18-gcccore-11.3.0-bare
#module load java
module load bedtools/2.30.0-gcc-10.3.0
#module load bamtools
#bsmap requires samtools < 1.0.0
# module load samtools/0.1.18
PATH=~/software/bsmap-2.74/samtools:$PATH

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAY_INDEX} eg 2 from {list}
ID="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p BSMAPratio
mkdir -p BSMAPratio_genome_mC_single_average
mkdir -p BSMAPratio_genome_mC_average_of_average
mkdir -p TempOut
#mkdir -p OnTargetCoverage
mkdir -p ConversionRate

########## Run #################
        # The required input is all in the folder bsmapped_filtered from the prior step
        ########################
        # extract methylation information using bsmap tool methratio.py
        # python ~/software/bsmap-2.74/methratio.py \
        # -o BSMAPratio/${ID}_methratio.txt \
        # -d ${genome_reference} \
        # -u \
        # -z \
        # -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

# output
#chr     pos     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower        CI_upper
#Chr01   74      -       AAGCC   0.000   7.00    0       7       0       0       0.000   0.354
#Chr01   158     -       TAGCC   0.000   1.00    0       1       0       0       0.000   0.793
#Chr01   174     -       CTGAG   0.000   2.00    0       2       0       0       0.000   0.658
#Chr01   176     -       GAGTT   0.000   2.00    0       2       0       0       0.000   0.658
#Chr01   186     -       AAGAT   0.000   4.00    0       4       0       0       0.000   0.490
#Chr01   195     +       ATCTA   0.000   1.00    0       1       3       3       0.000   0.793
#Chr01   201     +       TTCAG   0.000   1.00    0       1       3       3       0.000   0.793

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
        # CTG, CGG, CAG
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
                      else if(($3=="-" && $4~/^CAG../ ) || ($3=="+" &&  $4~/^..CTG/))
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

        # Output (no header)
        # adds end and makes it zero-based bed
        # parses nt context to get CG, CHG or CHH
        # keep the rest of the columns
        #chr     start  end     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower        CI_upper

        #Chr01   73      74      -       CHH     0.000   7.00    0       7       0       0       0.000   0.354
        #Chr01   157     158     -       CHH     0.000   1.00    0       1       0       0       0.000   0.793
        #Chr01   173     174     -       CHG     0.000   2.00    0       2       0       0       0.000   0.658
        #Chr01   175     176     -       CHH     0.000   2.00    0       2       0       0       0.000   0.658
        #Chr01   185     186     -       CHH     0.000   4.00    0       4       0       0       0.000   0.490
        #Chr01   194     195     +       CHH     0.000   1.00    0       1       3       3       0.000   0.793
        #Chr01   200     201     +       CHG     0.000   1.00    0       1       3       3       0.000   0.793


        if [ "$make_subcontext" == "yes" ]
        then
        awk -F$"\\t" "$awk_make_subcontext_bed" \
        "BSMAPratio/${ID}_methratio.txt" > "BSMAPratio/${ID}_BSMAP_out_subcontext.txt"
        fi
        ########################
        #For genome browser

        # begGraph ratio files for tdfs
        awk_make_bedGraph='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $8/$9*100, $5, $8, $9
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

        # pipe bedGraph to split by context (use dash to read from sdtin)
        # per context
        awk -F$"\\t" "$awk_make_bedGraph" \
        "BSMAPratio/${ID}_BSMAP_out.txt" | \
        awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" -

        # Summarise genome-wide methylation by suming C and CT per context for whole genome
        # This will reduce bias from low coverage sites (that are either 0 or 1) in low coverage datasets
        # use awk array per context ($5)
        awk -F$"\\t" -v ID=$ID 'BEGIN {OFS = FS} (NR>1){
          C_context[$5]+=$8; CT_context[$5]+=$9; next} END {
          for (i in C_context) print ID, i, C_context[i], CT_context[i], C_context[i]/CT_context[i]*100 > "BSMAPratio_genome_mC_single_average/"ID"_BSMAP_out_summary.txt" }
          ' "BSMAPratio/${ID}_BSMAP_out.txt"

        # Summarise genome-wide methylation by suming C and effCT per context for whole genome
        # This will reduce bias from low coverage sites (that are either 0 or 1) in low coverage datasets
        # use awk array per context ($5)
        awk -F$"\\t" -v ID=$ID 'BEGIN {OFS = FS} (NR>1){
          C_context[$5]+=$8; CT_context[$5]+=$7; next} END {
          for (i in C_context) print ID, i, C_context[i], CT_context[i], C_context[i]/CT_context[i]*100 > "BSMAPratio_genome_mC_single_average/"ID"_BSMAP_out_summary_eff-CT.txt" }
          ' "BSMAPratio/${ID}_BSMAP_out.txt"

        # add a step as above but taking the average of the average per cytosine ie average of ratio column
        # need to think though best awk method!

        if [ "$make_subcontext" == "yes" ]
        then
        # per sub-context
        awk -F$"\\t" "$awk_make_bedGraph" \
        "BSMAPratio/${ID}_BSMAP_out_subcontext.txt" | \
        awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_subcontext" -


        # per sub-context genome-wide methylation
        awk -F$"\\t" -v ID=$ID 'BEGIN {OFS = FS} (NR>1){
          C_context[$5]+=$8; CT_context[$5]+=$9; next} END {
          for (i in C_context) print ID, i, C_context[i], CT_context[i], C_context[i]/CT_context[i]*100 > "BSMAPratio_genome_mC_single_average/"ID"_BSMAP_out_subcontext_summary.txt" }
          ' "BSMAPratio/${ID}_BSMAP_out_subcontext.txt"

        fi

        #Make bigWigs per context
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CG.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_CG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CHG.bedGraph" ${chrom_sizes_file} \
        "BSMAPratio/${ID}_BSMAP_out_CHG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CHH.bedGraph" ${chrom_sizes_file} \
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
        if [ "$keep_bedgraph" == "no" ]
        then
        rm -rv BSMAPratio/${ID}*.bedGraph
        fi
        

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

                # Output (no header)
                # ratio calculated using C_count/CT_count
                #chr     start  end   ratio  context strand
                #Chr01   157     158     0       CHH     -
                #Chr01   173     174     0       CHG     -
                #Chr01   175     176     0       CHH     -
                #Chr01   185     186     0       CHH     -
                #Chr01   194     195     0       CHH     +
                #Chr01   200     201     0       CHG     +


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

                ## Make bigWigs per context
                # bigWigs for each individual induction with its own genome reference and plasmid
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_+_CG.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_+_CG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_+_CHG.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_+_CHG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_+_CHH.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_+_CHH.bigWig"

                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_-_CG.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_-_CG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_-_CHG.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_-_CHG.bigWig"
                bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_-_CHH.bedGraph" ${chrom_sizes_file} \
                "BSMAPratio/${ID}_BSMAP_out_-_CHH.bigWig"

        #remove bedGraph it is large and not really required
        # keep bigWigs?
        if [ "$keep_bedgraph" == "no" ]
        then
        rm -rv BSMAPratio/${ID}*.bedGraph
        fi

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
