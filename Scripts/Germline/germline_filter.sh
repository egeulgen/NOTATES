#!/bin/bash

##################################################
## Project: NOTATES
## Script purpose: Filter germline SNP/Indels
## Date: Oct 27, 2019
## Author: Ege Ulgen
##################################################


cd ./Germline
echo "##################################### Selecting Germline SNPs    " $(date)
$GATK SelectVariants -R $genome -V ./raw.snps.indels.vcf.gz \
	 --select-type-to-include SNP -O ./raw_snps.vcf.gz

echo "##################################### Filtering Germline SNPs    " $(date)
$GATK VariantFiltration -R $genome -V ./raw_snps.vcf.gz \
	--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filter-name "HQ_SNP_filter" -O ./germline_snps.vcf.gz

rm ./raw_snps.vcf
rm ./raw_snps.vcf.idx

echo "################################## Selecting Germline nonSNPs    " $(date)
### using --select-type-to-exclude SNP, in order to filter MIXED type with InDels
$GATK SelectVariants -R $genome -V ./raw.snps.indels.vcf.gz \
	--select-type-to-include INDEL -O ./raw_indels.vcf.gz

echo "################################## Filtering Germline nonSNPs    " $(date)
$GATK VariantFiltration -R $genome -V ./raw_indels.vcf.gz \
	--filter-expression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0" \
	--filter-name "HQ_InDel_filter" -O ./germline_indels.vcf.gz

rm ./raw_indels.vcf
rm ./raw_indels.vcf.idx

echo "########### Merging Germline SNP and InDels into a single VCF    " $(date)
$JAVA "$GATK3" -T CombineVariants -R $genome --variant ./germline_snps.vcf.gz \
	--variant ./germline_indels.vcf.gz -o ./filtered_germline_variants.vcf.gz \
	--genotypemergeoption UNSORTED

rm ./germline_snps.vcf
rm ./germline_snps.vcf.idx
rm ./germline_indels.vcf
rm ./germline_indels.vcf.idx

cd ..