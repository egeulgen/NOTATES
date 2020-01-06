#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Script for mapping (bwa),  
## Cleaning (picard), SAM>BAM (picard), fixing 
## mate info.(picard) and  marking duplicates 
## per lane. Combination of all anes (if >1 lane) 
## and BQSR (GATK)
## Date: Dec 6, 2019
## Author: Ege Ulgen
##################################################

analysisID=$1
sample_name=$2
sample_type=$3
kit_name=$4

cd ./$sample_name
lanes=$(cat ./lanes.txt)

####### Alignment, Clean SAM, SAM to BAM conversion, Fix Mate Information, and 
####### Mark Duplicates - by lane
mkdir QC

for lane in $lanes
do
	## Read Group Information
	readGroup="@RG\\tID:"$analysisID"\tPL:ILLUMINA\tLB:"$kit_name"\tSM:"$sample_name""

	## Alignment
	echo '############'"$sample_name"': Aligning reads from lane: '$lane "    " $(date)
	bwa mem -M -t "$num_threads" -R "$readGroup" "$genome" \
		"$lane"_R1.fastq.gz "$lane"_R2.fastq.gz > "$lane".sam

	rm "$lane"_R1.fastq.gz
	rm "$lane"_R2.fastq.gz

	## SAM to BAM
	echo '##########'"$sample_name"': SAM to BAM & Sort for lane: '$lane "    " $(date)
	$JAVA $PICARD SortSam SORT_ORDER=coordinate \
		INPUT="$lane".sam OUTPUT="$lane".bam CREATE_INDEX=true \
		USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
	
	rm "$lane".sam

	## Fix Mate Information
	echo '####'"$sample_name"': Fixing mate information for lane: '$lane "    " $(date)
	$JAVA $PICARD FixMateInformation SO=coordinate \
		INPUT="$lane".bam OUTPUT="$lane".fixed.bam ADD_MATE_CIGAR=TRUE \
		USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
	
	rm "$lane".bam "$lane".bai

	## Mark Duplicates - first pass by lane
	echo '##########'"$sample_name"': Marking duplicates of lane: '$lane "    " $(date)
	$JAVA $PICARD MarkDuplicates INPUT="$lane".fixed.bam \
		OUTPUT="$lane".marked.bam \
		METRICS_FILE=QC/"$lane"_MarkDup_metrics.txt CREATE_INDEX=true \
		USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
	
	rm "$lane".fixed.bam
done

####### Mark Duplicates - combine all lanes
if [ $(wc -w <<< "$lanes") != 1 ]
	then
	echo '######'"$sample_name"': Marking duplicates and merging BAM files    ' $(date)

	perlane=(${lanes// / })
	bams=("${perlane[@]/%/.marked.bam}")
	bais=("${perlane[@]/%/.marked.bai}")

	input=("${bams[@]/#/INPUT=}")

	$JAVA $PICARD MarkDuplicates "${input[@]}" OUTPUT="$sample_type".marked.bam \
		METRICS_FILE=QC/MarkDup_metrics.txt CREATE_INDEX=true \
		USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
	
	rm ${bams[@]} ${bais[@]}
else
	echo '#########'"$sample_name"': Renaming the marked BAM file (only a single lane)'
	
	mv "$lanes".marked.bam "$sample_type".marked.bam
	mv "$lanes".marked.bai "$sample_type".marked.bai
	mv QC/"$lanes"_MarkDup_metrics.txt QC/MarkDup_metrics.txt
fi

################################### Quality score recalibration
echo '##############################################'"$sample_name"': BQSR    ' $(date)
$GATK BaseRecalibrator -R "$genome" \
	--intervals "$Target_Intervals" --interval-padding 100 \
	-I "$sample_type".marked.bam \
	--known-sites "$dbSNP" --known-sites "$Mills_1kG" --known-sites "$ThousandG" \
	-O '../'"$sample_type"'.recal_data.table'

$GATK ApplyBQSR -R "$genome" \
	-I "$sample_type".marked.bam \
	--bqsr-recal-file ../"$sample_type".recal_data.table \
	-O ../"$sample_type".final.bam

rm ../"$sample_type".recal_data.table "$sample_type".marked.bam "$sample_type".marked.bai

cd ..