#!/bin/bash
################################################################################
######################### NeuroOncology Technologies ###########################
###################### Whole-Exome Sequencing Pipeline #########################
########################### Ege Ulgen, Dec 2019 ################################
################################################################################
#set -ueo pipefail

patientID=$1
normal_name=$2
tumor_name=$3
kit_name=$4
tumor_type=$5
primary_cond=$6
tumor_sample=$7

### Set main_dir to the directory where this script is located
main_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

### Source the configuration file
source "$main_dir"/notates.config "verbose" "$kit_name"

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
echo 
echo "The chosen patient ID is: ""$patientID"
echo "The germline sample ID is: ""$normal_name"
echo "The tumor sample ID  is: ""$tumor_name"
echo 
echo "The tumor type is: ""$tumor_type"
echo "The tumor sample type is: ""$tumor_sample"

################################################################################
############################## Merging fastq files #############################
################################################################################
echo "Merging fastq files (if necessary)"
####### Merge Seperate Files for normal
bash "$scripts_dir"/merge_fastqs.sh "$normal_name"

####### Merge Seperate Files for Tumor
bash "$scripts_dir"/merge_fastqs.sh "$tumor_name"

################################################################################
######################## Examine Sequence Read Quality #########################
################################################################################
echo "############################################## Running FASTQC    " $(date)
####### Merge Seperate Files for normal
bash "$scripts_dir"/fastqc.sh "$normal_name"

####### Merge Seperate Files for Tumor
bash "$scripts_dir"/fastqc.sh "$tumor_name"

################################################################################
######################### Mapping and Pre-processing ########################### 
################################################################################

####### Mapping and preprocessing for normal
bash "$scripts_dir"/mapping_preprocessing.sh "$patientID" "$normal_name" "normal" "$kit_name"

####### Mapping and preprocessing for tumor
bash "$scripts_dir"/mapping_preprocessing.sh "$patientID" "$tumor_name" "tumor" "$kit_name"

################################################################################
############################## Variant Calling #################################
################################################################################

########################################## HC ##################################
echo "############## Running Haplotype Caller for Germline Variants    " $(date)
mkdir Germline
$GATK HaplotypeCaller \
	-R "$genome" \
	-I normal.final.bam \
	--dbsnp "$dbSNP" \
	--intervals "$Target_Intervals" --interval-padding 100 \
	-O Germline/raw.snps.indels.vcf.gz

################################### Mutect #####################################
echo "############################################# Running MuTect2    " $(date)
mkdir Mutect2
$GATK Mutect2 \
	-R "$genome" \
	-I tumor.final.bam \
	-I normal.final.bam -normal "$normal_name" \
	--germline-resource "$gnomad_vcf" \
	--intervals "$Target_Intervals" --interval-padding 100 \
	-O Mutect2/Mutect_raw.vcf.gz

################################################################################
############################ Variant Filtering #################################
################################################################################

##################### Germline Variant Filtering ###############################
bash "$scripts_dir"/Germline/germline_filter.sh

##################### Somatic Variant Filtering ###############################
$GATK GetPileupSummaries \
	-R "$genome" \
	-I tumor.final.bam \
	-V "$small_exac_common" \
	--intervals "$small_exac_common" \
	-O tumor_pileups.table

$GATK GetPileupSummaries \
	-R "$genome" \
	-I normal.final.bam \
	-V "$small_exac_common" \
	--intervals "$small_exac_common" \
	-O normal_pileups.table

$GATK CalculateContamination \
	-I tumor_pileups.table \
	-matched normal_pileups.table \
	--tumor-segmentation segments.tsv \
	-O contamination.table

$GATK FilterMutectCalls \
	-R "$genome" \
	-V Mutect2/Mutect_raw.vcf.gz \
	--contamination-table contamination.table \
	--tumor-segmentation segments.tsv \
	-O Mutect2/Mutect_filt_once.vcf.gz

$GATK CollectSequencingArtifactMetrics \
	-I tumor.final.bam \
	-O artifact_metrics \
	-R $genome

## G/T: OxoG
## C/T: FFPE
$GATK FilterByOrientationBias \
	-AM 'G/T' \
	-AM 'C/T' \
	-V Mutect2/Mutect_filt_once.vcf.gz \
	-P artifact_metrics.pre_adapter_detail_metrics \
	-O Mutect2/Mutect_filt_twice.vcf.gz

$GATK SelectVariants -R "$genome" \
	-V Mutect2/Mutect_filt_twice.vcf.gz \
	--exclude-filtered \
	-O Mutect2/Filtered_mutect.vcf.gz

################################################################################
############################ Variant Annoation #################################
################################################################################
source "$oncotator_activate"
mkdir Oncotator

echo "################################# Annotating Somatic variants    " $(date)
oncotator -v --input_format=VCF --db-dir "$oncotator_ds" --tx-mode CANONICAL \
	-c "$oncotator_ds"/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
	Mutect2/Filtered_mutect.vcf.gz \
	Oncotator/annotated.sSNVs.tsv hg19

echo "################################ Annotating Germline Variants    " $(date)
oncotator -v --input_format=VCF --db-dir "$oncotator_ds" --tx-mode CANONICAL \
	-c "$oncotator_ds"/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
	Germline/filtered_germline_variants.vcf.gz \
	Oncotator/annotated.germline_SNVs.tsv hg19

deactivate

################################################################################
################################# ExomeCNV #####################################
################################################################################
bash "$scripts_dir"/ExomeCNV/ExomeCNV.sh $normal_name $tumor_name

################################################################################
################################## THetA #######################################
################################################################################
mkdir -p THetA/output

echo "################################## Preparing input for THetA     " $(date)
## ExomeCNV output to THetA input
bash "$THetA"/bin/CreateExomeInput -s ExomeCNV/CNV.segment.copynumber.txt \
	-t tumor.final.bam -n normal.final.bam \
	--FA $genome --EXON_FILE "$Target_Intervals" --QUALITY 30 --DIR THetA
rm tumor.final.pileup normal.final.pileup

echo "############################################### Running THetA    " $(date)
bash "$THetA"/bin/RunTHetA THetA/CNV.input \
	--DIR THetA/output --NUM_PROCESSES "$num_threads"
# bash "$THetA"/bin/RunTHetA THetA/CNV.input --TUMOR_FILE THetA/tumor_SNP.txt \
# 	--NORMAL_FILE THetA/normal_SNP.txt --DIR THetA/output --NUM_PROCESSES "$num_threads"

################################################################################
#################################### QC ########################################
################################################################################
## Alignment Summary Metrics
$JAVA $PICARD CollectAlignmentSummaryMetrics \
	METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$genome \
	INPUT=normal.final.bam OUTPUT="$normal_name"/QC/alignment_summary_metrics.txt

$JAVA $PICARD CollectAlignmentSummaryMetrics \
	METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$genome \
	INPUT=tumor.final.bam OUTPUT="$tumor_name"/QC/alignment_summary_metrics.txt

## QC Wrapper
Rscript "$scripts_dir"/QC_table_prep.R "$normal_name" "$tumor_name"

################################################################################
################################# NOTATES ######################################
################################################################################
echo "######################## Running R script for Germline Report    " $(date)
Rscript "$scripts_dir"/Germline/Germline.R

echo "######################## Running R script for NOTATES v4.1       " $(date)
Rscript "$scripts_dir"/NOTATES/run_NOTATES.R $tumor_type

echo "######################## Running R script for Pathway Enrichment " $(date)
Rscript "$scripts_dir"/pathway_enrichment.R

echo "######################## Running R script for DeConstructSigs    " $(date)
Rscript "$scripts_dir"/DeConstructSigs.R $scripts_dir $patientID

echo "######################## Running R script for MSIpred            " $(date)
Rscript "$scripts_dir"/MSIpred_prep.R $patientID
python "$scripts_dir"/MSIpred_analysis.py $simple_repeats $exome_length

echo "######################## Creating Report                         " $(date)
Rscript "$scripts_dir"/create_report.R $patientID $scripts_dir \
	$exome_length $Target_Intervals $tumor_type $primary_cond $tumor_sample

echo "######################## Finished 							   " $(date)
exit 0