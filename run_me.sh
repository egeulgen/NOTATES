#!/bin/bash
################################################################################
######################### NeuroOncology Technologies ###########################
###################### Whole-Exome Sequencing Pipeline #########################
########################## Ege Ulgen, March 2019 ###############################
################################################################################
#set -ueo pipefail

patientID=$1
normal_name=$2
tumor_name=$3
cap_kit=$4
tumor_type=$5

### Set main_dir to the directory where this script is located
main_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

### Source the configuration file
source "$main_dir"/configurations.cfg "verbose" $cap_kit

################################################################################
############################## Checkpoint for backup ###########################
################################################################################

read -p "Have you backed up original fastq files?(Y/N): "

if [ "$REPLY" != "Y" ]
	then 
	printf "Back up original files before proceeding\n"
	exit 1
fi

################################################################################
############################ Check Patient ID ##################################
######################## Normal ID and Tumor ID ################################
################################################################################

path_to_patient=$(find . -type d -name $patientID -print)
if [ ${#path_to_patient} == 0 ]
	then 
	printf "FOLDER NOT FOUND! :\n  Make sure that you are in the correct 
	directory\nand the patient ID is correct\n"
	exit 2
fi
cd $path_to_patient

if [ ! -d "$normal_name" ]
	then
	printf "FOLDER NOT FOUND! :\n  Make sure that the normal ID is correct\n"
	exit 3
fi

if [ ! -d "$tumor_name" ]
	then
	printf "FOLDER NOT FOUND! :\n  Make sure that the tumor ID is correct\n"
	exit 3
fi

#### Display the chosen patient ID and Normal ID, Tumor ID
echo ""
echo "The chosen patient ID is: ""$patientID"
echo "The germline sample ID is: ""$normal_name"
echo "The tumor sample ID  is: ""$tumor_name"
echo ""
echo "The tumor type is: ""$tumor_type"

################################################################################
############################## Merging fastq files #############################
################################################################################
echo "Merging fastq files (if necessary)"
####### Merge Seperate Files for normal
bash "$scripts_dir"/merge_fastqs.sh $normal_name

####### Merge Seperate Files for Tumor
bash "$scripts_dir"/merge_fastqs.sh $tumor_name

################################################################################
######################## Examine Sequence Read Quality #########################
################################################################################
echo "############################################## Running FASTQC    " $(date)
####### Merge Seperate Files for normal
bash "$scripts_dir"/fastqc.sh $normal_name

####### Merge Seperate Files for Tumor
bash "$scripts_dir"/fastqc.sh $tumor_name

################################################################################
######################### Mapping and Pre-processing ########################### 
################################################################################

####### Mapping and preprocessing for normal
bash "$scripts_dir"/mapping_preprocessing.sh $normal_name "normal"

####### Mapping and preprocessing for tumor
bash "$scripts_dir"/mapping_preprocessing.sh $tumor_name "tumor"

################################################################################
############################## Variant Calling #################################
################################################################################

########################################## HC ##################################
echo "############## Running Haplotype Caller for Germline Variants    " $(date)
mkdir ./Germline/
$GATK HaplotypeCaller -R $genome -I normal.final.bam --dbsnp $dbSNP \
	--intervals $Bait_Intervals --interval-padding 100 \
	-O ./Germline/raw.snps.indels.vcf

################################### Mutect #####################################
echo "############################################# Running MuTect2    " $(date)
mkdir ./Mutect
$GATK Mutect2 -R $genome \
	-I tumor.final.bam -tumor $tumor_name -I normal.final.bam -normal $normal_name \
	--germline-resource $gnomad_vcf --af-of-alleles-not-in-resource 0.0000025 \
	--intervals $Bait_Intervals --interval-padding 100 \
	-O ./Mutect/Mutect_raw.vcf.gz -bamout tumor_normal_m2.bam 

################################################################################
############################ Variant Filtering #################################
################################################################################

##################### Germline Variant Filtering ###############################
bash "$scripts_dir"/Germline/germline_filter.sh

##################### Somatic Variant Filtering ###############################
$GATK GetPileupSummaries \
-I tumor.final.bam \
-V $small_exac_common \
-O tumor_getpileupsummaries.table

$GATK GetPileupSummaries \
-I normal.final.bam \
-V $small_exac_common \
-O normal_getpileupsummaries.table

$GATK CalculateContamination \
-I tumor_getpileupsummaries.table -matched normal_getpileupsummaries.table \
-O tumor_calculatecontamination.table

$GATK FilterMutectCalls \
-V ./Mutect/Mutect_raw.vcf.gz \
--contamination-table tumor_calculatecontamination.table \
-O ./Mutect/Mutect_filt_once.vcf.gz

$GATK CollectSequencingArtifactMetrics \
-I tumor.final.bam \
-O artifact_metrics \
-R $genome

$GATK FilterByOrientationBias \
-AM 'G/T' \
-AM 'C/T' \
-V ./Mutect/Mutect_filt_once.vcf.gz \
-P artifact_metrics.pre_adapter_detail_metrics \
-O ./Mutect/Mutect_filt_twice.vcf.gz

$GATK SelectVariants -R $genome \
-V ./Mutect/Mutect_filt_twice.vcf.gz \
--exclude-filtered \
-O ./Mutect/Filtered_mutect.vcf.gz

################################################################################
############################ Variant Annoation #################################
################################################################################
source "$oncotator_activate"
mkdir ./Oncotator

echo "################################# Annotating Somatic variants    " $(date)
oncotator -v --input_format=VCF --db-dir "$oncotator_ds" --tx-mode CANONICAL \
	-c "$oncotator_ds"/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
	./Mutect/Filtered_mutect.vcf.gz ./Oncotator/annotated.sSNVs.tsv hg19

echo "################################ Annotating Germline Variants    " $(date)
oncotator -v --input_format=VCF --db-dir "$oncotator_ds" --tx-mode CANONICAL \
	-c "$oncotator_ds"/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
	./Germline/filtered_germline_variants.vcf \
	./Oncotator/annotated.germline_SNVs.tsv hg19

deactivate

################################################################################
########################### Germline Report ####################################
################################################################################
echo "######################## Running R script for Germline Report    " $(date)
Rscript "$scripts_dir"/Germline/Germline.R

################################################################################
################################# ExomeCNV #####################################
################################################################################
bash "$scripts_dir"/ExomeCNV/ExomeCNV.sh $normal_name $tumor_name

################################################################################
################################## THetA #######################################
################################################################################
mkdir -p ./THetA/output

echo "################################## Preparing input for THetA     " $(date)
## ExomeCNV output to THetA input
bash "$THetA"/bin/CreateExomeInput -s ./ExomeCNV/CNV.segment.copynumber.txt \
	-t tumor.final.bam -n normal.final.bam \
	--FA $genome --EXON_FILE $Bait_Intervals --QUALITY 30 --DIR ./THetA
rm tumor.final.pileup normal.final.pileup

echo "############################################### Running THetA    " $(date)
# bash "$THetA"/bin/RunTHetA THetA/CNV.input \
# 	--DIR ./THetA/output --NUM_PROCESSES 8
bash "$THetA"/bin/RunTHetA THetA/CNV.input --TUMOR_FILE THetA/tumor_SNP.txt \
	--NORMAL_FILE THetA/normal_SNP.txt --DIR ./THetA/output --NUM_PROCESSES 8

################################################################################
#################################### QC ########################################
################################################################################
## Insert-size Metrics
$JAVA $PICARD CollectInsertSizeMetrics \
	HISTOGRAM_FILE=$normal_name/QC/insert_size_histogram.pdf \
	INPUT=normal.final.bam OUTPUT=$normal_name/QC/insert_size_metrics.txt

$JAVA $PICARD CollectInsertSizeMetrics \
	HISTOGRAM_FILE=$tumor_name/QC/insert_size_histogram.pdf \
	INPUT=tumor.final.bam OUTPUT=$tumor_name/QC/insert_size_metrics.txt

## Depth of coverage over target intervals
$JAVA "$resources_dir""/Tools/GenomeAnalysisTK.jar" -T DepthOfCoverage -R $genome -I normal.final.bam \
	-o $normal_name/QC/primary_target_coverage \
	--intervals $Bait_Intervals --interval_padding 100 \
	-ct 1 -ct 5 -ct 10 -ct 25 -ct 50 -ct 100

$JAVA "$resources_dir""/Tools/GenomeAnalysisTK.jar" -T DepthOfCoverage -R $genome -I tumor.final.bam \
	-o $tumor_name/QC/primary_target_coverage \
	--intervals $Bait_Intervals --interval_padding 100 \
	-ct 1 -ct 5 -ct 10 -ct 25 -ct 50 -ct 100


## Alignment Summary Metrics
$JAVA $PICARD CollectAlignmentSummaryMetrics \
	METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$genome \
	INPUT=normal.final.bam OUTPUT=$normal_name/QC/alignment_summary_metrics.txt

$JAVA $PICARD CollectAlignmentSummaryMetrics \
	METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$genome \
	INPUT=tumor.final.bam OUTPUT=$tumor_name/QC/alignment_summary_metrics.txt

## FlagStat
samtools flagstat normal.final.bam > $normal_name/QC/flagstat_metrics.txt
samtools flagstat tumor.final.bam > $tumor_name/QC/flagstat_metrics.txt

## QC Wrapper
Rscript "$scripts_dir"/QC.R $normal_name $tumor_name

################################################################################
################################# NOTATES ######################################
################################################################################

echo "######################## Running R script for NOTATES v4.1       " $(date)
Rscript "$scripts_dir"/NOTATESv4.1/run_NOTATES.R $tumor_type

echo "######################## Running R script for Pathway Enrichment " $(date)
Rscript "$scripts_dir"/pathway_enrichment.R

echo "######################## Running R script for DeConstructSigs    " $(date)
Rscript "$scripts_dir"/DeConstructSigs.R $scripts_dir $patientID

echo "######################## Running R script for MSIpred            " $(date)
Rscript "$scripts_dir"/MSIpred_prep.R $patientID
python "$scripts_dir"/MSIpred_analysis.py $data_sources_dir $exome_length

echo "######################## Creating Report					       " $(date)
Rscript "$scripts_dir"/create_report.R $patientID $scripts_dir $exome_length $tumor_type

echo "######################## Finished 							   " $(date)
exit 0