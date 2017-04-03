#!/bin/bash
name=$1

cd ./"$name"
## Assuming Illumina fastq naming format
## <sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits)>.fastq.gz
	
export lanes=$(ls | grep 'L[0-9]*' -o | uniq)
echo "$lanes" >> lanes.txt

for lane in $lanes
do
	echo "Merging fastq files for: "$name" lane: "$lane"    " $(date)
	ls | grep $lane"_R1_"|xargs cat >> $lane'_R1.fastq.gz'
	ls | grep $lane"_R2_"|xargs cat >> $lane'_R2.fastq.gz'

	ls | grep $lane"_R1_"|xargs rm
	ls | grep $lane"_R2_"|xargs rm
done

cd ..