#!/bin/bash
name=$1

cd ./"$name"

lanes=$(ls | grep 'L[0-9]*' -o | uniq)
echo "$lanes" >> lanes.txt

for i in $lanes
do
	echo "Merging fastq files for: "$name" lane: "$i
	ls | grep $i"_R1_"|xargs cat >> $i'_R1.fastq.gz'
	ls | grep $i"_R2_"|xargs cat >> $i'_R2.fastq.gz'

	ls | grep $i"_R1_"|xargs rm
	ls | grep $i"_R2_"|xargs rm
done

cd ..

exit 0