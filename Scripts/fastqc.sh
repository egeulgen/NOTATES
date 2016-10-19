#!/bin/bash
name=$1

echo "FASTQC for ""$name    " $(date)
cd ./"$name"
mkdir QC

lanes=$(cat ./lanes.txt )

for i in $lanes
do
	echo "$i"
	fastqc --nogroup $i'_R1.fastq.gz' $i'_R2.fastq.gz' --outdir=./QC
done
cd ..

exit 0