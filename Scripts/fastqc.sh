#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Wrapper for FASTQC analysis per lane
## Date: Dec 6, 2019
## Author: Ege Ulgen
##################################################

name=$1

echo "################################## Running FASTQC for ""$name    " $(date)
cd ./"$name"
mkdir QC

lanes=$(cat ./lanes.txt)

for lane in $lanes
do
	echo "$lane"
	fastqc --nogroup "$lane"_R1.fastq.gz "$lane"_R2.fastq.gz --outdir=./QC
done
cd ..