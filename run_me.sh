#!/bin/bash
################################################################################
######################### NeuroOncology Technologies ###########################
###################### Whole-Exome Sequencing Pipeline #########################
########################## Ege Ulgen, July 2016 ################################
################################################################################

### Set main_dir to the directory where this script is located
main_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$main_dir"/configurations.cfg "verbose"

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
######################### Normal ID and Tumor ID################################
################################################################################
patientID=$1
if [ ${#patientID} == 0 ]
	then
	read -p "Enter the patient ID: " patientID
fi

path_to_patient=$(find . -type d -name $patientID -print)
if [ ${#path_to_patient} == 0 ]
	then 
	printf "FOLDER NOT FOUND! :\n  Make sure that the patient ID is correct\n"
	exit 2
fi
cd $path_to_patient

normal_name=$2
if [ ${#normal_name} == 0 ]
	then
	read -p "Enter the normal ID: " normal_name
fi
if [ ! -d "$normal_name" ]
	then
	printf "FOLDER NOT FOUND! :\n  Make sure that the normal ID is correct\n"
	exit 3
fi

tumor_name=$3
if [ ${#tumor_name} == 0 ]
	then
	read -p "Enter the tumor ID: " tumor_name
fi
if [ ! -d "$tumor_name" ]
	then
	printf "FOLDER NOT FOUND! :\n  Make sure that the tumor ID is correct\n"
	exit 3
fi

#### Display the chosen patient ID and Normal ID, Tumor ID
echo "The chosen patient ID is: ""$patientID"
echo "The germline sample ID is: ""$normal_name"
echo "The tumor sample ID  is: ""$tumor_name"
echo ""

################################################################################
############################## Merging fastq files #############################
################################################################################
echo "Merging fastq files"
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
########################## Contamination Estimation ############################ 
################################################################################
echo "##################### Conpair: Creating pileup file for tumor    " $(date)
${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py -B tumor.final.bam \
	-O tumor.conpair.pileup --reference $genome \
	--gatk $GATK --markers $conpair_markers_bed

echo "#################### Conpair: Creating pileup file for normal    " $(date)
${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py -B normal.final.bam \
	-O normal.conpair.pileup --reference $genome \
	--gatk $GATK --markers $conpair_markers_bed


echo "################################# Conpair: Verify concordance    " $(date)
${CONPAIR_DIR}/scripts/verify_concordance.py \
	-T tumor.conpair.pileup -N normal.conpair.pileup \
	--outfile conpair_concordance.txt --markers $conpair_markers_txt


echo "############################# Conpair: Estimate contamination    " $(date)
${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py \
	-T tumor.conpair.pileup -N normal.conpair.pileup \
	--outfile conpair_contam.txt --markers $conpair_markers_txt

rm tumor.conpair.pileup normal.conpair.pileup

Rscript "$scripts_dir"'/conpair_parser.R' $normal_name $tumor_name

################################################################################
############################## Variant Calling #################################
################################################################################

########################################## HC ##################################
echo "############## Running Haplotype Caller for Germline Variants    " $(date)
mkdir ./Germline/
$JAVA $GATK -T HaplotypeCaller -R $genome -I normal.final.bam --dbsnp $dbSNP \
	-stand_call_conf 30 -stand_emit_conf 10 \
	--intervals $Capture --interval_padding 100 \
	-o ./Germline/raw.snps.indels.vcf -nct 8

################################### MuTect #####################################
echo "############################################# Running MuTect2    " $(date)
mkdir ./MuTect
$JAVA $GATK -T MuTect2 -R $genome \
	-I:tumor tumor.final.bam -I:normal normal.final.bam \
	--dbsnp $dbSNP --cosmic $COSMIC \
	-contaminationFile contamination.txt \
	--intervals $Capture --interval_padding 100 \
	-o ./MuTect/MuTect_calls.vcf -nct 8

################################################################################
############################ Variant Filtering #################################
################################################################################

##################### Germline Variant Filtering ###############################
bash "$scripts_dir"/Germline/germline_filter.sh

################################################################################
############################ Variant Annoation #################################
################################################################################
source "$oncotator_activate"
mkdir ./Oncotator

echo "################################# Annotating Somatic variants    " $(date)
oncotator -v --input_format=VCF --db-dir "$oncotator_ds" --tx-mode CANONICAL \
	-c "$oncotator_ds"/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
	./MuTect/mutect_calls.vcf ./Oncotator/annotated.sSNVs.tsv hg19

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
bash "$scripts_dir"/ExomeCNV/ExomeCNV.sh $normal_name

################################################################################
################################## THetA #######################################
################################################################################
mkdir -p ./THetA/output

echo "################################## Preparing input for THetA     " $(date)
## ExomeCNV output to THetA input
bash "$THetA"/bin/CreateExomeInput -s ./ExomeCNV/CNV.segment.copynumber.txt \
	-t tumor.final.bam -n normal.final.bam \
	--FA $genome --EXON_FILE $Capture --QUALITY 30 --DIR ./THetA

echo "############################################### Running THetA    " $(date)
bash "$THetA"/bin/RunTHetA THetA/CNV.input --TUMOR_FILE THetA/tumor_SNP.txt \
	--NORMAL_FILE THetA/normal_SNP.txt --DIR ./THetA/output --NUM_PROCESSES 8

exit 0