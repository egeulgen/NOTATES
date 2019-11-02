#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Prepare necesarry data and 
## identify SCNAs using ExomeCNV
## Date: Oct 27, 2019
## Author: Ege Ulgen
##################################################

normal_name=$1
tumor_name=$2

####### Coverage files
echo "########################## Prep for ExomeCNV - Coverage files     " $(date)
mkdir -p ./ExomeCNV/DepthOfCoverage
echo "############################### Normal: Creating coverage file    " $(date)
$JAVA "$resources_dir""/Tools/GenomeAnalysisTK.jar" \
	-T DepthOfCoverage \
	-omitBaseOutput -omitLocusTable \
	-R $genome \
	-I normal.final.bam \
	--intervals $Bait_Intervals \
	-ct 1 -ct 5 -ct 10 -ct 25 -ct 50 -ct 100 \
	-o ./ExomeCNV/DepthOfCoverage/normal.coverage

echo "############################### Tumor: Creating coverage file     " $(date)
$JAVA "$resources_dir""/Tools/GenomeAnalysisTK.jar" \
	-T DepthOfCoverage \
	-omitBaseOutput -omitLocusTable \
	-R $genome \
	-I tumor.final.bam \
	--intervals $Bait_Intervals \
	-ct 1 -ct 5 -ct 10 -ct 25 -ct 50 -ct 100 \
	-o ./ExomeCNV/DepthOfCoverage/tumor.coverage

####### BAF files
echo "############################### Prep for ExomeCNV - BAF files    " $(date)
mkdir -p ./ExomeCNV/baf

## Get HQ-filtered Het SNP sites in Normal
expr1='vc.getGenotype("'"$normal_name"'").isHet()'
expr2='vc.getGenotype("'"$normal_name"'").getPhredScaledQual()>=30.0'
expr3='vc.getGenotype("'"$normal_name"'").getDP()>=20'

$GATK SelectVariants -R $genome \
	-V ./Germline/filtered_germline_variants.vcf -select-type SNP \
	--selectExpressions "$expr1" \
	--selectExpressions "$expr2" \
	--selectExpressions "$expr3" \
	--exclude-filtered --restrict-alleles-to BIALLELIC \
	-O ./ExomeCNV/baf/normal_HQ_SNPs.vcf

## Call variants at the same HQ-filtered Het SNP sites in Tumor
# HC
$GATK HaplotypeCaller -R $genome -I tumor.final.bam \
	--intervals ./ExomeCNV/baf/normal_HQ_SNPs.vcf \
	-O ./ExomeCNV/baf/raw_tumor_HC.vcf

# Filter
$GATK VariantFiltration -R $genome -V ./ExomeCNV/baf/raw_tumor_HC.vcf \
	--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filter-name "HQ_SNP_filter" -O ./ExomeCNV/baf/filtered_tumor_HC.vcf

rm ./ExomeCNV/baf/raw_tumor_HC.vcf ./ExomeCNV/baf/raw_tumor_HC.vcf.idx

# Select variants
expr1='vc.getGenotype("'"$tumor_name"'").getPhredScaledQual() >= 30.0'
expr2='vc.getGenotype("'"$tumor_name"'").getDP() >= 20'

$GATK SelectVariants -R $genome \
	-V ./ExomeCNV/baf/filtered_tumor_HC.vcf -select-type SNP \
	--selectExpressions "$expr1" \
	--selectExpressions "$expr2" \
	--exclude-filtered --restrict-alleles-to BIALLELIC \
	-O ./ExomeCNV/baf/tumor_HQ_SNPs.vcf

rm ./ExomeCNV/baf/filtered_tumor_HC.vcf ./ExomeCNV/baf/filtered_tumor_HC.vcf.idx

Rscript "$scripts_dir"'/ExomeCNV/VCF_parser_BAF.R'

####### ExomeCNV
echo "############################################# Running R script for ExomeCNV    " $(date)
Rscript "$scripts_dir"'/ExomeCNV/ExomeCNV.R' $read_length

exit 0