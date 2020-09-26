#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Prepare necesarry data and 
## identify SCNAs using ExomeCNV
## Date: Sep 23, 2020
## Author: Ege Ulgen
##################################################

normal_name=$1
tumor_name=$2

####### Coverage files
echo "########################## Prep for ExomeCNV - Coverage files     " $(date)
mkdir -p ./ExomeCNV/DepthOfCoverage
echo "############################### Normal: Creating coverage file    " $(date)
gatk3 "$java_options"\
	-T DepthOfCoverage \
	-omitBaseOutput -omitLocusTable \
	-R "$genome" \
	-I normal.final.bam \
	--intervals "$Target_Intervals" \
	-ct 1 -ct 5 -ct 10 -ct 25 -ct 50 -ct 100 \
	-o ./ExomeCNV/DepthOfCoverage/normal.coverage

echo "############################### Tumor: Creating coverage file     " $(date)
gatk3 "$java_options"\
	-T DepthOfCoverage \
	-omitBaseOutput -omitLocusTable \
	-R "$genome" \
	-I tumor.final.bam \
	--intervals "$Target_Intervals" \
	-ct 1 -ct 5 -ct 10 -ct 25 -ct 50 -ct 100 \
	-o ./ExomeCNV/DepthOfCoverage/tumor.coverage

####### BAF files
echo "############################### Prep for ExomeCNV - BAF files    " $(date)
mkdir -p ./ExomeCNV/baf

## Get HQ-filtered Het SNP sites in Normal
expr1='vc.getGenotype("'"$normal_name"'").isHet()'
expr2='vc.getGenotype("'"$normal_name"'").getPhredScaledQual()>=30.0'
expr3='vc.getGenotype("'"$normal_name"'").getDP()>=20'

gatk SelectVariants \
	--java-options "$java_options" \
	-R "$genome" \
	-V ./Germline/filtered_germline_variants.vcf.gz \
	-select-type SNP \
	--selectExpressions "$expr1" \
	--selectExpressions "$expr2" \
	--selectExpressions "$expr3" \
	--exclude-filtered \
	--restrict-alleles-to BIALLELIC \
	-O ./ExomeCNV/baf/normal_HQ_SNPs.vcf

## Call variants at the same HQ-filtered Het SNP sites in Tumor
# HC
gatk HaplotypeCaller \
	--java-options "$java_options" \
	-R "$genome" \
	-I tumor.final.bam \
	--intervals ./ExomeCNV/baf/normal_HQ_SNPs.vcf \
	-O ./ExomeCNV/baf/raw_tumor_HC.vcf.gz

# Filter
gatk VariantFiltration \
	--java-options "$java_options" \
	-R $genome \
	-V ./ExomeCNV/baf/raw_tumor_HC.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O ./ExomeCNV/baf/filtered_tumor_HC.vcf.gz

rm ./ExomeCNV/baf/raw_tumor_HC.vcf.gz*

# Select variants
expr1='vc.getGenotype("'"$tumor_name"'").getPhredScaledQual() >= 30.0'
expr2='vc.getGenotype("'"$tumor_name"'").getDP() >= 20'

gatk SelectVariants \
	--java-options "$java_options" \
	-R $genome \
	-V ./ExomeCNV/baf/filtered_tumor_HC.vcf.gz \
	-select-type SNP \
	--selectExpressions "$expr1" \
	--selectExpressions "$expr2" \
	--exclude-filtered \
	--restrict-alleles-to BIALLELIC \
	-O ./ExomeCNV/baf/tumor_HQ_SNPs.vcf

rm ./ExomeCNV/baf/filtered_tumor_HC.vcf.gz*

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_R

Rscript "$scripts_dir"/ExomeCNV/VCF_parser_BAF.R

####### ExomeCNV
echo "############################################# Running R script for ExomeCNV    " $(date)
Rscript "$scripts_dir"/ExomeCNV/ExomeCNV.R "$read_length" "$num_threads"

exit 0