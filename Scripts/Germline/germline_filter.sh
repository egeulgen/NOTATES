#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Filter germline SNP/Indels
## Date: Sep 18, 2020
## Author: Ege Ulgen
##################################################


cd Germline
echo "##################################### Selecting Germline SNPs    " $(date)
$GATK SelectVariants \
	-R "$genome" \
	-V raw.snps.indels.vcf.gz \
	--select-type-to-include SNP \
	--select-type-to-exclude SYMBOLIC \
	--restrict-alleles-to BIALLELIC \
	--exclude-non-variants \
	-O raw_snps.vcf.gz

echo "##################################### Filtering Germline SNPs    " $(date)
$GATK VariantFiltration \
	-R "$genome" \
	-V raw_snps.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O germline_snps.vcf.gz

rm raw_snps.vcf.gz
rm raw_snps.vcf.gz.tbi

echo "################################## Selecting Germline nonSNPs    " $(date)
$GATK SelectVariants \
	-R "$genome" \
	-V raw.snps.indels.vcf.gz \
	--select-type-to-include INDEL \
	--select-type-to-exclude SYMBOLIC \
	--restrict-alleles-to BIALLELIC \
	-O raw_indels.vcf.gz

echo "################################## Filtering Germline nonSNPs    " $(date)
$GATK VariantFiltration \
	-R "$genome" \
	-V raw_indels.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 200.0" --filter-name "FS200" \
	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	-O germline_indels.vcf.gz

rm raw_indels.vcf.gz
rm raw_indels.vcf.gz.tbi

echo "########### Merging Germline SNP and InDels into a single VCF    " $(date)
$JAVA "$GATK3" \
	-T CombineVariants \
	-R "$genome" \
	--variant germline_snps.vcf.gz \
	--variant germline_indels.vcf.gz \
	-o combined_germline_variants.vcf.gz \
	--genotypemergeoption UNSORTED

rm germline_snps.vcf.gz* germline_indels.vcf.gz*

$GATK SelectVariants -R "$genome" \
	-V combined_germline_variants.vcf.gz \
	--exclude-filtered \
	-O filtered_germline_variants.vcf.gz

cd ..
