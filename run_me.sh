#!/bin/bash
################################################################################
######################### NeuroOncology Technologies ###########################
###################### Whole-Exome Sequencing Pipeline #########################
########################### Ege Ulgen, Sep 2020 ################################
################################################################################
# set -ueo pipefail

export CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

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

conda activate NOTATES_main
################################################################################
######################## Examine Sequence Read Quality #########################
################################################################################
echo "############################################## Running FASTQC    " $(date)
####### FASTQC for normal
bash "$scripts_dir"/fastqc.sh "$normal_name"

####### FASTQC for Tumor
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
gatk HaplotypeCaller \
	--java-options "$java_options" \
	-R "$genome" \
	-I normal.final.bam \
	--dbsnp "$dbSNP" \
	--intervals "$Target_Intervals" --interval-padding 100 \
	-O Germline/raw.snps.indels.vcf.gz

################################### Mutect #####################################
echo "############################################# Running MuTect2    " $(date)
mkdir Mutect2
gatk Mutect2 \
	--java-options "$java_options" \
	-R "$genome" \
	-I tumor.final.bam \
	-I normal.final.bam -normal "$normal_name" \
	--germline-resource "$gnomad_vcf" \
	--intervals "$Target_Intervals" --interval-padding 100 \
	--f1r2-tar-gz Mutect2/f1r2.tar.gz \
	-O Mutect2/Mutect_raw.vcf.gz

################################################################################
############################ Variant Filtering #################################
################################################################################

##################### Germline Variant Filtering ###############################
bash "$scripts_dir"/Germline/germline_filter.sh

##################### Somatic Variant Filtering ###############################
bash "$scripts_dir"/somatic_filter.sh

################################################################################
############################ Variant Annoation #################################
################################################################################
mkdir Funcotator

echo "################################ Annotating Somatic Variants    " $(date)
gatk Funcotator \
	--java-options "$java_options" \
	-R "$genome" \
	-V Mutect2/Filtered_mutect.vcf.gz \
	--ref-version "$ref_version" \
	--data-sources-path "$funcotator_ds_somatic" \
	--output Funcotator/annotated_somatic.maf \
	--output-file-format MAF \
	--annotation-default Center:NOT \
	--annotation-default tumor_barcode:"$tumor_name" \
	--annotation-default normal_name:"$normal_name" \
	--transcript-list "$funcotator_ds_somatic"/transcriptList.exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt

echo "################################ Annotating Germline Variants   " $(date)
gatk Funcotator \
	--java-options "$java_options" \
	-R "$genome" \
	-V Germline/filtered_germline_variants.vcf.gz \
	--ref-version "$ref_version" \
	--data-sources-path "$funcotator_ds_germline" \
	--output Funcotator/annotated_germline.maf \
	--output-file-format MAF \
	--annotation-default Center:NOT \
	--annotation-default normal_name:"$normal_name" 

################################################################################
################################# ExomeCNV #####################################
################################################################################
bash "$scripts_dir"/ExomeCNV/ExomeCNV.sh $normal_name $tumor_name

################################################################################
################################## THetA #######################################
################################################################################

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_python
mkdir -p THetA/output

echo "################################## Preparing input for THetA     " $(date)
## ExomeCNV output to THetA input
bash "$THetA"/bin/CreateExomeInput -s ExomeCNV/CNV.segment.copynumber.txt \
	-t tumor.final.bam -n normal.final.bam \
	--FA $genome --EXON_FILE "$Target_Intervals" --QUALITY 30 --DIR THetA
rm tumor.final.pileup normal.final.pileup

echo "############################################### Running THetA    " $(date)
bash "$THetA"/bin/RunTHetA THetA/CNV.input \
	--DIR THetA/output --NUM_PROCESSES "$num_threads" --MIN_FRAC 0.005
# bash "$THetA"/bin/RunTHetA THetA/CNV.input --TUMOR_FILE THetA/tumor_SNP.txt \
# 	--NORMAL_FILE THetA/normal_SNP.txt --DIR THetA/output --NUM_PROCESSES "$num_threads"
conda deactivate

################################################################################
################################## MSIpred #####################################
################################################################################
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_R
echo "######################## Running R script for MSIpred            " $(date)
Rscript "$scripts_dir"/MSIpred_prep.R $patientID
conda deactivate

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_python
python "$scripts_dir"/MSIpred_analysis.py $simple_repeats $exome_length
conda deactivate

conda deactivate

################################################################################
################################ Translocations ################################
################################################################################
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_main

bash "$scripts_dir"/SV_calling.sh

conda deactivate

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_R

Rscript "$scripts_dir"/process_DELLY.R

conda deactivate


################################################################################
#################################### QC ########################################
################################################################################
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_main
## Alignment Summary Metrics
picard CollectAlignmentSummaryMetrics \
	METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$genome \
	INPUT=normal.final.bam OUTPUT="$normal_name"/QC/alignment_summary_metrics.txt \
	USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

picard CollectAlignmentSummaryMetrics \
	METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$genome \
	INPUT=tumor.final.bam OUTPUT="$tumor_name"/QC/alignment_summary_metrics.txt \
	USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
conda deactivate

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate NOTATES_R
# QC Wrapper
Rscript "$scripts_dir"/QC_table_prep.R "$normal_name" "$tumor_name"

################################################################################
################################# NOTATES ######################################
################################################################################
echo "######################## Running R script for Germline Report    " $(date)
Rscript "$scripts_dir"/Germline/Germline.R

echo "######################## Running R script for NOTATES Report     " $(date)
Rscript "$scripts_dir"/NOTATES/run_NOTATES.R $tumor_type

echo "######################## Running R script for Pathway Enrichment " $(date)
Rscript "$scripts_dir"/pathway_enrichment.R

echo "######################## Running R script for DeConstructSigs    " $(date)
Rscript "$scripts_dir"/DeConstructSigs.R $scripts_dir $patientID

echo "######################## Creating Report                         " $(date)
Rscript "$scripts_dir"/create_report.R $patientID $scripts_dir \
	$exome_length $Target_Intervals $tumor_type $primary_cond $tumor_sample

conda deactivate
echo "######################## Finished 							   " $(date)
exit 0