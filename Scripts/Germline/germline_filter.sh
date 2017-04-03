#!/bin/bash

cd ./Germline
echo "##################################### Selecting Germline SNPs    " $(date)
$JAVA $GATK -T SelectVariants -R $genome -V ./raw.snps.indels.vcf \
	-selectType SNP -o ./raw_snps.vcf

echo "##################################### Filtering Germline SNPs    " $(date)
$JAVA $GATK -T VariantFiltration -R $genome -V ./raw_snps.vcf \
	--filterExpression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "HQ_SNP_filter" -o ./germline_snps.vcf

rm ./raw_snps.vcf
rm ./raw_snps.vcf.idx

echo "############################# Selecting Germline non-SNP variants" $(date)
### using --selectTypeToExclude SNP, in order to filter MIXED type with InDels
$JAVA $GATK -T SelectVariants -R $genome -V ./raw.snps.indels.vcf \
	--selectTypeToExclude SNP -o ./raw_indels.vcf

echo "############################# Filtering Germline non-SNP variants" $(date)
$JAVA $GATK -T VariantFiltration -R $genome -V ./raw_indels.vcf \
	--filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0" \
	--filterName "HQ_InDel_filter" -o ./germline_indels.vcf

rm ./raw_indels.vcf
rm ./raw_indels.vcf.idx

echo "########### Merging Germline SNP and InDels into a single VCF    " $(date)
$JAVA $GATK -T CombineVariants -R $genome --variant ./germline_snps.vcf \
	--variant ./germline_indels.vcf -o ./filtered_germline_variants.vcf \
	--genotypemergeoption UNSORTED

rm ./germline_snps.vcf
rm ./germline_snps.vcf.idx
rm ./germline_indels.vcf
rm ./germline_indels.vcf.idx

cd ..