#!/bin/bash

normal_name=$1

####### Coverage files
echo "########################## Prep for ExomeCNV - Coverage files    " $(date)
mkdir -p ./ExomeCNV/DepthOfCoverage
echo "############################### Normal: Creating overage file    " $(date)
$JAVA $GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $genome \
	-I normal.final.bam --intervals $Capture \
	-o ./ExomeCNV/DepthOfCoverage/normal.coverage

echo "############################### Tumor: Creating coverage file    " $(date)
$JAVA $GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $genome \
	-I tumor.final.bam --intervals $Capture \
	-o ./ExomeCNV/DepthOfCoverage/tumor.coverage

####### BAF files
echo "############################### Prep for ExomeCNV - BAF files    " $(date)
mkdir -p ./ExomeCNV/baf

## Get HQ-filtered Het SNP sites
expr='vc.getGenotype("'"$normal_name"'").isHet()'
$JAVA $GATK -T SelectVariants -R $genome \
	-V ./Germline/filtered_germline_variants.vcf -selectType SNP \
	-select $expr -o ./ExomeCNV/baf/HQ_het_sites.vcf --excludeFiltered

skip=$(($(sed -n "/#CHROM/=" ./ExomeCNV/baf/HQ_het_sites.vcf) + 1))
cat ./ExomeCNV/baf/HQ_het_sites.vcf | tail -n +$skip | awk '{FS="\t";OFS="\t";print $1,$2,$2}' > ./ExomeCNV/baf/HQ_het_sites.bed
rm ./ExomeCNV/baf/HQ_het_sites.vcf ./ExomeCNV/baf/HQ_het_sites.vcf.idx

bam-readcount -q 20 -b 10 -f $genome -w 1 -l ./ExomeCNV/baf/HQ_het_sites.bed ./normal.final.bam > ./ExomeCNV/baf/normal_bamreadcount.txt
bam-readcount -q 20 -b 10 -f $genome -w 1 -l ./ExomeCNV/baf/HQ_het_sites.bed ./tumor.final.bam > ./ExomeCNV/baf/tumor_bamreadcount.txt

Rscript "$scripts_dir"'/ExomeCNV/BAM_readcount_to_BAF.R' 20

####### ExomeCNV
echo "############################################# Running R script for ExomeCNV    " $(date)
Rscript "$scripts_dir"'/ExomeCNV/ExomeCNV.R' $(pwd)

exit 0