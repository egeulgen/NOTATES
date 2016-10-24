#!/bin/bash

normal_name=$1
tumor_name=$2

####### Coverage files
echo "########################## Prep for ExomeCNV - Coverage files    " $(date)
mkdir -p ./ExomeCNV/DepthOfCoverage
echo "############################### Normal: Creating overage file    " $(date)
$JAVA $GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $genome \
	-I normal.final.bam --intervals $Bait_Intervals \
	-o ./ExomeCNV/DepthOfCoverage/normal.coverage

echo "############################### Tumor: Creating coverage file    " $(date)
$JAVA $GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $genome \
	-I tumor.final.bam --intervals $Bait_Intervals \
	-o ./ExomeCNV/DepthOfCoverage/tumor.coverage

####### BAF files
echo "############################### Prep for ExomeCNV - BAF files    " $(date)
mkdir -p ./ExomeCNV/baf

## Get HQ-filtered Het SNP sites in Normal
expr1='vc.getGenotype("'"$normal_name"'").isHet()'
expr2='vc.getGenotype("'"$normal_name"'").getPhredScaledQual() >= 30.0'
expr3='vc.getGenotype("'"$normal_name"'").getDP() >= 20'

$JAVA $GATK -T SelectVariants -R $genome \
	-V ./Germline/filtered_germline_variants.vcf -selectType SNP \
	--selectexpressions "$expr1" \
	--selectexpressions "$expr2" --selectexpressions "$expr3" \
	--excludeFiltered -restrictAllelesTo BIALLELIC \
	-o ./ExomeCNV/baf/normal_HQ_SNPs.vcf

## Get the same HQ-filtered Het SNP sites in Tumor
# HC
$JAVA $GATK -T HaplotypeCaller -R $genome -I tumor.final.bam \
	-stand_call_conf 30 -stand_emit_conf 10 \
	--intervals ./ExomeCNV/baf/normal_HQ_SNPs.vcf \
	-o ./ExomeCNV/baf/raw_tumor_HC.vcf -nct 8

# Filter
$JAVA $GATK -T VariantFiltration -R $genome -V ./ExomeCNV/baf/raw_tumor_HC.vcf \
	--filterExpression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "HQ_SNP_filter" -o ./ExomeCNV/baf/filtered_tumor_HC.vcf

rm ./ExomeCNV/baf/raw_tumor_HC.vcf ./ExomeCNV/baf/raw_tumor_HC.vcf.idx

# Select variants
expr='vc.getGenotype("'"$tumor_name"'").getDP() >= 20'

$JAVA $GATK -T SelectVariants -R $genome \
	-V ./ExomeCNV/baf/filtered_tumor_HC.vcf -selectType SNP \
	--selectexpressions "$expr" \
	--excludeFiltered -restrictAllelesTo BIALLELIC \
	-o ./ExomeCNV/baf/tumor_HQ_SNPs.vcf

rm ./ExomeCNV/baf/filtered_tumor_HC.vcf ./ExomeCNV/baf/filtered_tumor_HC.vcf.idx

Rscript "$scripts_dir"'/ExomeCNV/VCF_parser_BAF.R'

####### ExomeCNV
echo "############################################# Running R script for ExomeCNV    " $(date)
Rscript "$scripts_dir"'/ExomeCNV/ExomeCNV.R'

exit 0