#!/bin/bash

name=$1
sample=$2

cd ./$name
lanes=$(cat ./lanes.txt)

####### Alignment, Clean SAM, SAM to BAM conversion, Fix Mate Information, and 
####### Mark Duplicates - by lane

flowcell=$(sed '2q;d' SampleSheet.csv | cut -f1 -d",")
for lane in $lanes
do
	## Read Group Information
	RG='@RG'
	RG=$RG'\tID:'$flowcell'.'$lane
	RG=$RG'\tSM:'$name
	RG=$RG'\tPL:ILLUMINA\tLB:library'

	## Alignment
	echo '############'"$name"': Aligning reads from lane: '$lane "    " $(date)
	bwa mem -M -t 16 -R $RG $genome \
		$lane'_R1.fastq.gz' $lane'_R2.fastq.gz' > $lane'.sam'

	rm $lane'_R1.fastq.gz'
	rm $lane'_R2.fastq.gz'

	## Clean SAM
	echo '################'"$name"': Cleaning SAM of lane: '$lane "    " $(date)
	$JAVA $PICARD CleanSam INPUT=$lane'.sam' OUTPUT=$lane'.clean.sam'
	
	rm $lane'.sam'

	## SAM to BAM
	echo '##########'"$name"': SAM to BAM & Sort for lane: '$lane "    " $(date)
	$JAVA $PICARD SortSam SORT_ORDER=coordinate \
		INPUT=$lane'.clean.sam' OUTPUT=$lane'.bam' CREATE_INDEX=true
	
	rm $lane'.clean.sam'

	## Fix Mate Information
	echo '####'"$name"': Fixing mate information for lane: '$lane "    " $(date)
	$JAVA $PICARD FixMateInformation SO=coordinate \
		INPUT=$lane'.bam' OUTPUT=$lane'.fixed.bam' ADD_MATE_CIGAR=TRUE
	
	rm $lane'.bam' $lane'.bai'

	## Mark Duplicates - first pass by lane
	echo '##########'"$name"': Marking duplicates of lane: '$lane "    " $(date)
	$JAVA $PICARD MarkDuplicates INPUT=$lane'.fixed.bam' \
		OUTPUT=$lane'.marked.bam' \
		METRICS_FILE='./QC/'$lane'_MarkDup_metrics.txt' CREATE_INDEX=true
	
	rm $lane'.fixed.bam'
done

####### Mark Duplicates - combine all lanes
if [ $(wc -w <<< "$lanes") != 1 ]
	then
	echo '######'"$name"': Marking duplicates and merging BAM files    ' $(date)

	perlane=(${lanes// / })
	bams=("${perlane[@]/%/.marked.bam}")
	bais=("${perlane[@]/%/.marked.bai}")

	input=("${bams[@]/#/INPUT=}")

	$JAVA $PICARD MarkDuplicates ${input[@]} OUTPUT="$sample"'.marked.bam' \
		METRICS_FILE='./QC/MarkDup_metrics.txt' CREATE_INDEX=true
	
	rm ${bams[@]} ${bais[@]}
else
	echo '#########'"$name"': Renaming the marked BAM file (only a single lane)'
	
	mv $lanes'.marked.bam' "$sample"'.marked.bam'
	mv $lanes'.marked.bai' "$sample"'.marked.bai'
	mv './QC/'"$lanes"'_MarkDup_metrics.txt' './QC/MarkDup_metrics.txt'
fi
rm lanes.txt

################################### Quality score recalibration
echo '##############################################'"$name"': BQSR    ' $(date)
$JAVA $GATK -T BaseRecalibrator -R $genome \
	--intervals $Bait_Intervals --interval_padding 100 \
	-I "$sample"'.marked.bam' \
	-knownSites $dbSNP -knownSites $Mills_1kG -knownSites $ThousandG \
	-o '../'"$sample"'.recal_data.table'  -nct 8 
$JAVA $GATK -T PrintReads -R $genome -I "$sample"'.marked.bam' \
	-BQSR '../'"$sample"'.recal_data.table' \
	-o '../'"$sample"'.final.bam' -nct 8
cd ..