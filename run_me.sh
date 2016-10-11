#!/bin/bash
##########################################################################################
############################## NeuroOncology Technologies ################################
########################### Whole-Exome Sequencing Pipeline ##############################
############################### Ege Ulgen, July 2016 #####################################
##########################################################################################

### Set main_dir to the directory where this script is located
main_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$main_dir"/configurations.cfg verbose

##########################################################################################
################################### Checkpoint for backup ################################ 
##########################################################################################

echo "Have you backed up original fastq files?(Y/N)"
read backup

if [ "$backup" != "Y" ]
	then 
	printf "Back up original files before proceeding\n"
	exit 1
fi

##########################################################################################
################################## Check Patient ID ######################################
##########################################################################################
echo "Enter the patient ID:"
read patientID

path=$(find . -type d -name $patientID -print)
if [ ${#path} == 0 ]
	then 
	printf "FOLDER NOT FOUND! :\n  Make sure that the patient ID is correct\n"
	exit 1
fi

cd $path
echo "Enter the normal ID:"
read normal_name

if [ ! -d "$normal_name" ]
	then
	printf "FOLDER NOT FOUND! :\n  Make sure that the ID is correct\n"
	exit 1
fi

echo "Enter the tumor ID:"
read tumor_name

if [ ! -d "$tumor_name" ]
	then
	printf "FOLDER NOT FOUND! :\n  Make sure that the ID is correct\n"
	exit 1
fi

#### Display the chosen patient ID and Normal ID, Tumor ID
echo "The chosen patient ID is: ""$patientID"
echo "The germline sample ID is: ""$normal_name"
echo "The tumor sample ID  is: ""$tumor_name"
echo ""
##########################################################################################
################################### Merging fastq files ##################################
########################################################################################## 
echo "Merging fastq files"
####### Merge Seperate Files for normal
bash "$scripts_dir"/merge_fastqs.sh $normal_name

####### Merge Seperate Files for Tumor
bash "$scripts_dir"/merge_fastqs.sh $tumor_name

##########################################################################################
############################# Examine Sequence Read Quality ##############################
##########################################################################################
echo "FASTQC"
####### Merge Seperate Files for normal
bash "$scripts_dir"/fastqc.sh $normal_name

####### Merge Seperate Files for Tumor
bash "$scripts_dir"/fastqc.sh $tumor_name

##########################################################################################
############################## Mapping and Pre-processing ################################ 
##########################################################################################

####### Mapping and preprocessing for normal
bash "$scripts_dir"/mapping_preprocessing.sh $main_dir $normal_name "normal"

####### Mapping and preprocessing for tumor
bash "$scripts_dir"/mapping_preprocessing.sh $main_dir $tumor_name "tumor"

##########################################################################################
################################### Variant Calling ######################################
##########################################################################################

########################################## HC ############################################
echo "############################################# Running Haplotype Caller for Germline Variants"
mkdir ./Germline/
$JAVA $GATK -T HaplotypeCaller -R $genome -I normal.final.bam --dbsnp $dbSNP -stand_call_conf 30 -stand_emit_conf 10 -L $Capture -o ./Germline/raw.snps.indels.vcf -nct 8
clear

######################################## MuTect ##########################################
echo "############################################# Running MuTect2"
mkdir ./MuTect
$JAVA $GATK -T MuTect2 -R $genome -I:tumor tumor.final.bam -I:normal normal.final.bam --dbsnp $dbSNP --cosmic $COSMIC -L $Capture -o ./MuTect/mutect_calls.vcf -nct 8
clear

##########################################################################################
################################# Variant Filtering ######################################
##########################################################################################

########################## Germline Variant Filtering ####################################
bash "$scripts_dir"/Germline/germline_filter.sh $main_dir

##########################################################################################
################################# Variant Annoation ######################################
##########################################################################################
source "$oncotator_activate"
mkdir ./Oncotator
Data_sources
echo "############################################# Annotating Somatic variants"
oncotator -v --input_format=VCF --db-dir "$oncotator_ds" --tx-mode CANONICAL -c "$oncotator_ds"/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt ./MuTect/mutect_calls.vcf ./Oncotator/annotated.sSNVs.tsv hg19

echo "############################################# Annotating Germline Variants"
oncotator -v --input_format=VCF --db-dir "$oncotator_ds" --tx-mode CANONICAL -c "$oncotator_ds"/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt ./Germline/filtered_germline_variants.vcf ./Oncotator/annotated.germline_SNVs.tsv hg19

deactivate

##########################################################################################
################################ Germline Report #########################################
##########################################################################################

echo "############################################# Running R script for Germline Report"
Rscript "$scripts_dir"/Germline/Germline.R $(pwd) $main_dir

##########################################################################################
###################################### ExomeCNV ##########################################
##########################################################################################

bash "$scripts_dir"/ExomeCNV/ExomeCNV.sh $main_dir $normal_name


exit 0