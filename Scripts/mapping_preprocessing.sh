#!/bin/bash

main_dir=$1
name=$2
sample=$3

source "$main_dir"/configurations.cfg 

cd ./$name
lanes=$(cat ./lanes.txt )

####### Alignment, Clean SAM, SAM to BAM conversion, Fix Mate Information, and Mark Duplicates - by lane
flowcell=$(sed '2q;d' SampleSheet.csv | cut -f1 -d",")
for i in $lanes
do
	## Alignment
	echo '#############################################'"$name"': Aligning reads from lane: '$i "    " $(date)
	bwa mem -M -t 16 -R '@RG\tID:'"$flowcell"'.'"$i"'\tSM:'"$name"'\tPL:ILLUMINA\tLB:library' $genome $i'_R1.fastq.gz' $i'_R2.fastq.gz' > $i'.sam'

	rm $i'_R1.fastq.gz'
	rm $i'_R2.fastq.gz'

	## Clean SAM
	echo '#############################################'"$name"': Cleaning SAM of lane: '$i "    " $(date)
	$JAVA $PICARD CleanSam INPUT=$i'.sam' OUTPUT=$i'.clean.sam'
	
	rm $i'.sam'

	## SAM to BAM
	echo '#############################################'"$name"': SAM to BAM & Sort of lane: '$i "    " $(date)
	$JAVA $PICARD SortSam SORT_ORDER=coordinate INPUT=$i'.clean.sam' OUTPUT=$i'.bam' CREATE_INDEX=true
	
	rm $i'.clean.sam'

	## Fix Mate Information
	echo '#############################################'"$name"': Fixing mate information of lane: '$i "    " $(date)
	$JAVA $PICARD FixMateInformation SO=coordinate INPUT=$i'.bam' OUTPUT=$i'.fixed.bam' ADD_MATE_CIGAR=TRUE
	
	rm $i'.bam' $i'.bai'

	## Mark Duplicates - first pass by lane
	echo '#############################################'"$name"': Marking duplicates of lane: '$i "    " $(date)
	$JAVA $PICARD MarkDuplicates INPUT=$i'.fixed.bam' OUTPUT=$i'.marked.bam' METRICS_FILE='./QC/'"$i"'_MarkDup_metrics.txt' CREATE_INDEX=true
	
	rm $i'.fixed.bam'
done

####### Mark Duplicates - combine all lanes
if [ $(wc -w <<< "$lanes") != 1 ]
	then
	echo '#############################################'"$name"': Marking duplicates and merging BAM files    ' $(date)

	perlane=(${lanes// / })
	bams=("${perlane[@]/%/.marked.bam}")
	bais=("${perlane[@]/%/.marked.bai}")

	input=("${bams[@]/#/INPUT=}")

	$JAVA $PICARD MarkDuplicates ${input[@]} OUTPUT="$sample"'.marked.bam' METRICS_FILE='./QC/MarkDup_metrics.txt' CREATE_INDEX=true
	
	rm ${bams[@]} ${bais[@]}
else
	echo '#############################################'"$name"': Renaming the marked BAM files, because only there is only a single lane'
	
	mv $lanes'.marked.bam' "$sample"'.marked.bam'
	mv $lanes'.marked.bai' "$sample"'.marked.bai'
	mv './QC/'"$lanes"'_MarkDup_metrics.txt' './QC/MarkDup_metrics.txt'
fi

################################### Quality score recalibration
echo '#############################################'"$name"': BQSR    ' $(date)
$JAVA $GATK -T BaseRecalibrator -R $genome -L $Capture -I "$sample"'.marked.bam' -knownSites $dbSNP -knownSites $Mills_1kG -knownSites $ThousandG -o '../'"$sample"'.recal_data.table'  -nct 8 
$JAVA $GATK -T PrintReads -R $genome -I "$sample"'.marked.bam' -BQSR '../'"$sample"'.recal_data.table'  -o '../'"$sample"'.final.bam' -nct 8

cd ..
exit 0