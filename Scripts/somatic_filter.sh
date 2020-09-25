#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Filter somatic SNV/Indels
## Date: Sep 24, 2020
## Author: Ege Ulgen
##################################################

gatk LearnReadOrientationModel \
	--java-options "$java_options" \
	-I Mutect2/f1r2.tar.gz \
	-O Mutect2/read-orientation-model.tar.gz

gatk GetPileupSummaries \
	--java-options "$java_options" \
	-R "$genome" \
	-I tumor.final.bam \
	-V "$small_exac_common" \
	--intervals "$small_exac_common" \
	-O tumor_pileups.table

gatk GetPileupSummaries \
	--java-options "$java_options" \
	-R "$genome" \
	-I normal.final.bam \
	-V "$small_exac_common" \
	--intervals "$small_exac_common" \
	-O normal_pileups.table

gatk CalculateContamination \
	--java-options "$java_options" \
	-I tumor_pileups.table \
	-matched normal_pileups.table \
	--tumor-segmentation segments.tsv \
	-O contamination.table

gatk FilterMutectCalls \
	--java-options "$java_options" \
	-R "$genome" \
	-V Mutect2/Mutect_raw.vcf.gz \
	--tumor-segmentation segments.tsv \
	--contamination-table contamination.table \
	--ob-priors Mutect2/read-orientation-model.tar.gz \
	-O Mutect2/Mutect_filt_applied.vcf.gz

gatk SelectVariants -R "$genome" \
	--java-options "$java_options" \
	-V Mutect2/Mutect_filt_applied.vcf.gz \
	--exclude-filtered \
	-O Mutect2/Filtered_mutect.vcf.gz
