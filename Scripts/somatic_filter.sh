#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Filter somatic SNV/Indels
## Date: Sep 18, 2020
## Author: Ege Ulgen
##################################################

$GATK LearnReadOrientationModel \
	-I Mutect2/f1r2.tar.gz \
	-O Mutect2/read-orientation-model.tar.gz

$GATK GetPileupSummaries \
	-R "$genome" \
	-I tumor.final.bam \
	-V "$small_exac_common" \
	--intervals "$small_exac_common" \
	-O tumor_pileups.table

$GATK GetPileupSummaries \
	-R "$genome" \
	-I normal.final.bam \
	-V "$small_exac_common" \
	--intervals "$small_exac_common" \
	-O normal_pileups.table

$GATK CalculateContamination \
	-I tumor_pileups.table \
	-matched normal_pileups.table \
	--tumor-segmentation segments.tsv \
	-O contamination.table

$GATK FilterMutectCalls \
	-R "$genome" \
	-V Mutect2/Mutect_raw.vcf.gz \
	--tumor-segmentation segments.tsv \
	--contamination-table contamination.table \
	--ob-priors Mutect2/read-orientation-model.tar.gz \
	-O Mutect2/Mutect_filt_applied.vcf.gz

$GATK SelectVariants -R "$genome" \
	-V Mutect2/Mutect_filt_applied.vcf.gz \
	--exclude-filtered \
	-O Mutect2/Filtered_mutect.vcf.gz