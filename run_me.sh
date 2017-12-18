#!/bin/bash
################################################################################
######################### NeuroOncology Technologies ###########################
###################### Whole-Exome Sequencing Pipeline #########################
########################### Ege Ulgen, March 2017 ##############################
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
######################## Normal ID and Tumor ID ################################
################################################################################
patientID=$1
if [ ${#patientID} == 0 ]
	then
	read -p "Enter the patient ID: " patientID
fi

path_to_patient=$(find . -type d -name $patientID -print)
if [ ${#path_to_patient} == 0 ]
	then 
	printf "FOLDER NOT FOUND! :\n  Make sure that you are in the correct 
	directory\nand the patient ID is correct\n"
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
echo ""
echo "The chosen patient ID is: ""$patientID"
echo "The germline sample ID is: ""$normal_name"
echo "The tumor sample ID  is: ""$tumor_name"
echo ""

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
	-stand_call_conf 30 \
	--intervals $Bait_Intervals --interval_padding 100 \
	-o ./Germline/raw.snps.indels.vcf -nct 8

################################### MuTect #####################################
echo "############################################# Running MuTect2    " $(date)
mkdir ./MuTect
$JAVA $GATK -T MuTect2 -R $genome \
	-I:tumor tumor.final.bam -I:normal normal.final.bam \
	--dbsnp $dbSNP --cosmic $COSMIC \
	-contaminationFile contamination.txt \
	--intervals $Bait_Intervals --interval_padding 100 \
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
$JAVA $GATK -T DepthOfCoverage -R $genome -I normal.final.bam \
	-o $normal_name/QC/primary_target_coverage \
	--intervals $Bait_Intervals --interval_padding 100 \
	-ct 1 -ct 5 -ct 10 -ct 25 -ct 50 -ct 100

$JAVA $GATK -T DepthOfCoverage -R $genome -I tumor.final.bam \
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

## Before and after plots for BQSR
$JAVA $GATK -T BaseRecalibrator -R $genome -I ./$normal_name/normal.marked.bam \
	-knownSites $dbSNP -knownSites $Mills_1kG -knownSites $ThousandG \
	-BQSR normal.recal_data.table -o normal.after_recal.table \
	--intervals $Bait_Intervals --interval_padding 100

$JAVA $GATK -T BaseRecalibrator -R $genome -I ./$tumor_name/tumor.marked.bam \
	-knownSites $dbSNP -knownSites $Mills_1kG -knownSites $ThousandG \
	-BQSR tumor.recal_data.table -o tumor.after_recal.table \
	--intervals $Bait_Intervals --interval_padding 100

rm ./$normal_name/normal.marked.bam 
rm ./$normal_name/normal.marked.bai
rm ./$tumor_name/tumor.marked.bam
rm ./$tumor_name/tumor.marked.bai

$JAVA $GATK -T AnalyzeCovariates -R $genome \
	-before normal.recal_data.table -after normal.after_recal.table \
	-plots ./$normal_name/BQSR_covariates.pdf

$JAVA $GATK -T AnalyzeCovariates -R $genome \
	-before tumor.recal_data.table -after tumor.after_recal.table \
	-plots ./$tumor_name/BQSR_covariates.pdf

## QC Wrapper
Rscript "$scripts_dir"/QC.R $normal_name $tumor_name

################################################################################
################################# NOTATES ######################################
################################################################################

echo "######################## Running R script for NOTATES v4.1       " $(date)
Rscript "$scripts_dir"/NOTATESv4.1/run_NOTATES.R

echo "######################## Running R script for Pathway Enrichment " $(date)
Rscript "$scripts_dir"/pathway_enrichment.R

echo "######################## Running R script for DeConstructSigs    " $(date)
Rscript "$scripts_dir"/DeConstructSigs.R $scripts_dir

echo "######################## Creating Report					       " $(date)
cp $scripts_dir/Report.Rmd ./Report.Rmd
Rscript "$scripts_dir"/create_report.R $patientID $scripts_dir
rm Report.Rmd
mv Report.pdf Report_"$patientID".pdf

echo "######################## Finished 							   " $(date)
exit 0