#!/bin/bash
name=$1

cd ./"$name"
## Assuming Illumina fastq naming format
## <sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits)>.fastq.gz
	
export lanes=$(ls | grep 'L[0-9]*' -o | uniq)
echo "$lanes" >> lanes.txt

for lane in $lanes
do
	export R1_files=$(ls | grep $lane"_R1_")
	if [ ${#R1_files[@]} != 1 ]
		then
		echo "Merging fastq files for: "$name" lane: "$lane" R1  " $(date)
		ls | grep $lane"_R1_"|xargs cat >> $lane'_R1.fastq.gz'
		ls | grep $lane"_R1_"|xargs rm
	else
		mv $R1_files $lane'_R1.fastq.gz'
	fi

	export R2_files=$(ls | grep $lane"_R2_")
	if [ ${#R2_files[@]} != 1 ]
		then
		echo "Merging fastq files for: "$name" lane: "$lane" R2  " $(date)
		ls | grep $lane"_R2_"|xargs cat >> $lane'_R2.fastq.gz'
		ls | grep $lane"_R2_"|xargs rm
	else
		mv $R2_files $lane'_R2.fastq.gz'
	fi
done

cd ..