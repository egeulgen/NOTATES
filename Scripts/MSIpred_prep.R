##################################################
## Project: NOTATES
## Script purpose: Prepare MAF file for MSIpred
## Analysis
## Date: Sep 22, 2020
## Author: Ege Ulgen
##################################################

args <- commandArgs(trailingOnly=TRUE)

dir.create("MSIpred")
funcotator_maf <- read.delim("Funcotator/annotated_somatic.maf", comment.char="#")
funcotator_maf$Tumor_Sample_Barcode <- args[1]
write.table(funcotator_maf, "MSIpred/somatic_maf.maf", sep = "\t", quote = FALSE, row.names = FALSE) 
