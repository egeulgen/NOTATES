#!/bin/bash

main_dir=$1
normal_name=$2

source "$main_dir"/configurations.cfg 

####### Coverage files
echo "############################################# Prep for ExomeCNV - Coverage files"
mkdir -p ./ExomeCNV/DepthOfCoverage
echo "############################ Normal: Creating coverage file"
$JAVA $GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $genome -I normal.final.bam -L $Capture -o ./ExomeCNV/DepthOfCoverage/normal.coverage

echo "############################ Tumor: Creating coverage file"
$JAVA $GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $genome -I tumor.final.bam -L $Capture -o ./ExomeCNV/DepthOfCoverage/tumor.coverage

####### BAF files
echo "############################################# Prep for ExomeCNV - BAF files"
mkdir -p ./ExomeCNV/baf

## Get HQ-filtered Het sites
expr='vc.getGenotype("'"$normal_name"'").isHet()'
$JAVA $GATK -T SelectVariants -R $genome -V ./Germline/filtered_germline_variants.vcf -selectType SNP -select $expr -o ./ExomeCNV/baf/HQ_het_sites.vcf --excludeFiltered

skip=$(($(sed -n "/#CHROM/=" ./ExomeCNV/baf/HQ_het_sites.vcf) + 1))
cat ./ExomeCNV/baf/HQ_het_sites.vcf | tail -n +$skip | awk '{FS="\t";OFS="\t";print $1,$2}' > ./ExomeCNV/baf/HQ_het_sites.bed
#cat ./ExomeCNV/baf/HQ_het_sites.vcf | tail -n +$skip | awk '{FS="\t";OFS="\t";print $1,$2-1,$2,$3}' > ./ExomeCNV/baf/HQ_het_sites.bed
rm ./ExomeCNV/baf/HQ_het_sites.vcf ./ExomeCNV/baf/HQ_het_sites.vcf.idx


## Normal BAF file
samtools mpileup -BQ0 -d10000000 -f $genome -l ./ExomeCNV/baf/HQ_het_sites.bed -o normal.pileup normal.final.bam
cat normal.pileup | "$main_dir"'/Scripts/ExomeCNV/mpileup2baf.pl' --min-reads 20 > ./ExomeCNV/baf/normal_baf.txt
rm normal.pileup

## Tumor BAF file
samtools mpileup -BQ0 -d10000000 -f $genome -l ./ExomeCNV/baf/HQ_het_sites.bed -o tumor.pileup tumor.final.bam
cat tumor.pileup | "$main_dir"'/Scripts/ExomeCNV/mpileup2baf.pl' --min-reads 20 > ./ExomeCNV/baf/tumor_baf.txt
rm tumor.pileup

####### ExomeCNV
echo "############################################# Running R script for ExomeCNV"
Rscript "$scripts_dir"'/ExomeCNV/ExomeCNV.R' $(pwd) $patientID --save

exit 0