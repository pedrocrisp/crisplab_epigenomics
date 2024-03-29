#!/bin/bash
#PBS -N 04-bigwigs
#PBS -r y
#PBS -m abej

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
module load samtools/1.10
module load bedtools/2.25.0

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAY_INDEX}p ${LIST})"

sample=$ID
sample_dir=analysis/${alignFolder}

outdir=analysis/bigWig_intermediate_${alignFolder}
mkdir -p ${outdir}

bigWigs_outdir=analysis/${alignFolder}_bigWigs
mkdir -p ${bigWigs_outdir}

# Condition statement: if library is stranded $3 == stranded and bigWigs are made for each strand.  If library is nonstranded $3 == nonstranded and nonstrandspecific bigWig is made.  If no strand info is specified, script will error
if [ "$strand" == "stranded_PE" ]
then

#split F and R reads into plus and minus strand taking into account PE
#http://seqanswers.com/forums/showthread.php?t=29399

#R1 forward strand
samtools view -f 99 -b $sample_dir/$sample.bam > ${outdir}/${sample}.R1F.bam

#R2 reverse strand
samtools view -f 147 -b $sample_dir/$sample.bam > ${outdir}/${sample}.R2R.bam

samtools merge -f ${outdir}/${sample}.forward.bam ${outdir}/${sample}.R1F.bam ${outdir}/${sample}.R2R.bam

#R1 reverse strand
samtools view -f 83 -b $sample_dir/$sample.bam > ${outdir}/${sample}.R1R.bam

#R2 forward strand
samtools view -f 163 -b ${sample_dir}/$sample.bam > ${outdir}/${sample}.R2F.bam

samtools merge -f ${outdir}/${sample}.reverse.bam ${outdir}/${sample}.R1R.bam ${outdir}/${sample}.R2F.bam

rm ${outdir}/${sample}*.R*.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bed

#make tdf - never use it... comment out
#echo "bedgraph to binary tiled data (.tdf) file"
#igvtools toTDF $outdir/${sample}.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
sort -k1,1 -k2,2n $outdir/${sample}.bedgraph > $outdir/${sample}.sorted.bedgraph
bedGraphToBigWig $outdir/${sample}.sorted.bedgraph $chrc_sizes $bigWigs_outdir/$sample.bw

sort -k1,1 -k2,2n $outdir/${sample}.plus.bg > $outdir/${sample}.sorted.plus.bg
bedGraphToBigWig $outdir/${sample}.sorted.plus.bg $chrc_sizes $bigWigs_outdir/$sample.plus.bw

sort -k1,1 -k2,2n $outdir/${sample}.minus.bg > $outdir/${sample}.sorted.minus.bg
bedGraphToBigWig $outdir/${sample}.sorted.minus.bg $chrc_sizes $bigWigs_outdir/$sample.minus.bw

elif [ "$strand" == "reverse_stranded_PE" ]
then


#split F and R reads into plus and minus strand taking into account PE
#http://seqanswers.com/forums/showthread.php?t=29399

#R1 forward strand
samtools view -f 99 -b $sample_dir/$sample.bam   > ${outdir}/${sample}.R1F.bam

#R2 reverse strand
samtools view -f 147 -b $sample_dir/$sample.bam   > ${outdir}/${sample}.R2R.bam

samtools merge -f ${outdir}/${sample}.forward.bam ${outdir}/${sample}.R1F.bam ${outdir}/${sample}.R2R.bam

#R1 reverse strand
samtools view -f 83 -b $sample_dir/$sample.bam   > ${outdir}/${sample}.R1R.bam

#R2 forward strand
samtools view -f 163 -b $sample_dir/$sample.bam   > ${outdir}/${sample}.R2F.bam

samtools merge -f ${outdir}/${sample}.reverse.bam ${outdir}/${sample}.R1R.bam ${outdir}/${sample}.R2F.bam

rm ${outdir}/${sample}*.R*.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bed

#make tdf - never use it... comment out
#echo "bedgraph to binary tiled data (.tdf) file"
#igvtools toTDF $outdir/${sample}.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
sort -k1,1 -k2,2n $outdir/${sample}.bedgraph > $outdir/${sample}.sorted.bedgraph
bedGraphToBigWig $outdir/${sample}.sorted.bedgraph $chrc_sizes $bigWigs_outdir/$sample.bw

sort -k1,1 -k2,2n $outdir/${sample}.plus.bg > $outdir/${sample}.sorted.plus.bg
bedGraphToBigWig $outdir/${sample}.sorted.plus.bg $chrc_sizes $bigWigs_outdir/$sample.plus.bw

sort -k1,1 -k2,2n $outdir/${sample}.minus.bg > $outdir/${sample}.sorted.minus.bg
bedGraphToBigWig $outdir/${sample}.sorted.minus.bg $chrc_sizes $bigWigs_outdir/$sample.minus.bw

elif [ "$strand" == "stranded_SE" ]
then

#R1 forward strand -f 0 does not work, all reads retained... "samtools view -f 0 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.forward.bam"
#not 100% sure we want everything that is NOT unmapped or reserse mapped, not sure what else coudl be included... of well
samtools view -F 4 -b $sample_dir/$sample.bam | samtools view -F 16 -b - > ${outdir}/${sample}.forward.bam

#R1 reverse strand
samtools view -f 16 -b $sample_dir/$sample.bam   > ${outdir}/${sample}.reverse.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph
# sort
bedSort $outdir/${sample}.bedgraph $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bg
# sort
bedSort $outdir/${sample}.minus.bg $outdir/${sample}.minus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bg
# sort
bedSort $outdir/${sample}.plus.bg $outdir/${sample}.plus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bed

#commented out becasue of issue on MSI, dont use these files anyway
#make tdf
#echo "bedgraph to binary tiled data (.tdf) file"
#igvtools toTDF $outdir/${sample}.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"

sort -k1,1 -k2,2n $outdir/${sample}.bedgraph > $outdir/${sample}.sorted.bedgraph
bedGraphToBigWig $outdir/${sample}.sorted.bedgraph $chrc_sizes $bigWigs_outdir/$sample.bw

sort -k1,1 -k2,2n $outdir/${sample}.plus.bg > $outdir/${sample}.sorted.plus.bg
bedGraphToBigWig $outdir/${sample}.sorted.plus.bg $chrc_sizes $bigWigs_outdir/$sample.plus.bw

sort -k1,1 -k2,2n $outdir/${sample}.minus.bg > $outdir/${sample}.sorted.minus.bg
bedGraphToBigWig $outdir/${sample}.sorted.minus.bg $chrc_sizes $bigWigs_outdir/$sample.minus.bw
###########

##

elif [ "$strand" == "stranded_SE_5prime" ]
then

#R1 forward strand -f 0 does not work, all reads retained... "samtools view -f 0 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.forward.bam"
#not 100% sure we want everything that is NOT unmapped or reserse mapped, not sure what else coudl be included... of well
samtools view -F 4 -b $sample_dir/$sample.bam | samtools view -F 16 -b - > ${outdir}/${sample}.forward.bam

#R1 reverse strand
samtools view -f 16 -b $sample_dir/$sample.bam   > ${outdir}/${sample}.reverse.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph
# sort
bedSort $outdir/${sample}.bedgraph $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bg
# sort
bedSort $outdir/${sample}.minus.bg $outdir/${sample}.minus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bg
# sort
bedSort $outdir/${sample}.plus.bg $outdir/${sample}.plus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus.bed

#commented out becasue of issue on MSI, dont use these files anyway
#make tdf
#echo "bedgraph to binary tiled data (.tdf) file"
#igvtools toTDF $outdir/${sample}.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
bedGraphToBigWig $outdir/${sample}.bedgraph $chrc_sizes $bigWigs_outdir/$sample.bigWig

bedGraphToBigWig $outdir/${sample}.plus.bg $chrc_sizes $bigWigs_outdir/$sample.plus.bigWig
bedGraphToBigWig $outdir/${sample}.minus.bg $chrc_sizes $bigWigs_outdir/$sample.minus.bigWig

### 5' analysis
### Note consider adding this to other sections, alternatively leaving it here will probably mean it only gets run for PARE data (single end stranded)
#stranded bedgraphs coverage of 5' position only
bedtools genomecov -5 -d -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus5.bed
#minus strand reads bedgraph
bedtools genomecov -5 -d -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus5.bed

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -5 -bga -scale -1 -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.minus5.bg
# sort
bedSort $outdir/${sample}.minus5.bg $outdir/${sample}.minus5.bg
#minus strand reads bedgraph
bedtools genomecov -5 -bga -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.plus5.bg
# sort
bedSort $outdir/${sample}.plus5.bg $outdir/${sample}.plus5.bg

bedGraphToBigWig $outdir/${sample}.plus5.bg $chrc_sizes $bigWigs_outdir/$sample.plus5.bigWig
bedGraphToBigWig $outdir/${sample}.minus5.bg $chrc_sizes $bigWigs_outdir/$sample.minus5.bigWig
###########

##

elif [ "$strand" == "reverse_stranded_SE" ]
then

#R1 forward strand -f 0 does not work, all reads retained... "samtools view -f 0 -b $sample_dir/$sample.bam   > $sample_dir/${sample}.forward.bam"
#not 100% sure we want everything that is NOT unmapped or reserse mapped, not sure what else coudl be included... of well
samtools view -F 4 -b $sample_dir/$sample.bam | samtools view -F 16 -b - > ${outdir}/${sample}.forward.bam

#R1 reverse strand
samtools view -f 16 -b $sample_dir/$sample.bam   > ${outdir}/${sample}.reverse.bam

####################

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -bga -split -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bg
#minus strand reads bedgraph
bedtools genomecov -bga -split -scale -1 -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bg

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
bedtools genomecov -d -split -ibam ${outdir}/${sample}.reverse.bam -g $chrc_sizes > $outdir/${sample}.plus.bed
#minus strand reads bedgraph
bedtools genomecov -d -split -scale -1 -ibam ${outdir}/${sample}.forward.bam -g $chrc_sizes > $outdir/${sample}.minus.bed

#make tdf - never use it... comment out
#echo "bedgraph to binary tiled data (.tdf) file"
#igvtools toTDF $outdir/${sample}.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
sort -k1,1 -k2,2n $outdir/${sample}.bedgraph > $outdir/${sample}.sorted.bedgraph
bedGraphToBigWig $outdir/${sample}.sorted.bedgraph $chrc_sizes $bigWigs_outdir/$sample.bw

sort -k1,1 -k2,2n $outdir/${sample}.plus.bg > $outdir/${sample}.sorted.plus.bg
bedGraphToBigWig $outdir/${sample}.sorted.plus.bg $chrc_sizes $bigWigs_outdir/$sample.plus.bw

sort -k1,1 -k2,2n $outdir/${sample}.minus.bg > $outdir/${sample}.sorted.minus.bg
bedGraphToBigWig $outdir/${sample}.sorted.minus.bg $chrc_sizes $bigWigs_outdir/$sample.minus.bw

##

elif [ "$strand" == "nonstranded" ]
then

echo "bam to bedgraph"
#non-stranded bedgraph
bedtools genomecov -bga -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bedgraph

#stranded bedgraphs with splicing and nt resolution - not using the '-strand +' flag because accounting for PE reads
#plus strand reads bedgraph
#could add this option to all steps above - output .bed file is about the same size as the bam ie liek 1.8 GB... quite big!
bedtools genomecov -d -split -ibam $sample_dir/$sample.bam -g $chrc_sizes > $outdir/${sample}.bed

#make tdf - never use it... comment out
#echo "bedgraph to binary tiled data (.tdf) file"
#igvtools toTDF $outdir/${sample}.bedgraph $outdir/$sample.tdf $chrc_sizes

#make bigWigs
echo "bam to bigWig"
sort -k1,1 -k2,2n $outdir/${sample}.bedgraph > $outdir/${sample}.sorted.bedgraph
bedGraphToBigWig $outdir/${sample}.sorted.bedgraph $chrc_sizes $bigWigs_outdir/$sample.bw

else
echo "ERROR: it has not been specificed whether library is stranded on not"

exit 1
fi

echo Done making beds
